# author: taewoo
# date: Apr 9 2025
# description: Fig. 2b-e - number of point driver mutations and tumor ploidy (mean copy number)
#   by diagnostic group (FTA, minimally invasive FTC, widely invasive FTC) and in FTC
#   non-invasive vs. invasive regions.
# Notes:
#   - All Wilcoxon rank-sum tests are exact. base R wilcox.test() falls back to the normal
#     approximation when ties are present (as with driver counts), so wx_exact() uses
#     coin::wilcox_test(distribution = "exact"), the conditional exact distribution with ties.
#   - The p-values printed on the figure (ANNOTATIONS) are checked against the computed values
#     with stopifnot(), so the figure and the test results cannot diverge.
#   - Segment files are read with read.table(..., comment.char = "") because their header line
#     starts with '#CHROM'.
#   - Points show per-patient means; tests use the first sample per patient (and per invasive status).

# SET DIRECTORY -----------------------------------------------------------

# Paths: run from the repository root; input data are read from data/ and figures are written to figures/.
CODE_DIR   <- normalizePath(".")
DATA_DIR   <- file.path(CODE_DIR, "data")
FIGUREPATH <- file.path(CODE_DIR, "figures")
dir.create(FIGUREPATH, showWarnings = FALSE)
setwd(FIGUREPATH)
# Library -----------------------------------------------------------------

suppressMessages(library(dplyr)); 
suppressMessages(library(readxl))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(grid))
suppressMessages(library(coin))   # exact tests with ties (see header)
suppressMessages(library(ggpubr))
library(dplyr)
library(purrr)
library(tidyr)

source(file.path(CODE_DIR, 'Thyroid_sources.R'))

# METADATA ----------------------------------------------------------------
metadata_file = file.path(DATA_DIR, "thyroid_sample_metadata.xlsx")
metadata <- read_sample_metadata(metadata_file)

# patient order (by Diagnosis + Ras)
patient_order <- metadata %>%
  group_by(Diagnosis, Ras, Patient) %>%
  summarise(max_Tip_total_mut = max(Tip_total_mut, na.rm = TRUE), .groups = "drop") %>%
  arrange(
    factor(Diagnosis, levels = c("FTA", "minimallyFTC", "widelyFTC")),
    factor(Ras, levels = c("F", "T")),
    max_Tip_total_mut
  ) %>%
  pull(Patient)

# process metadata and apply the final ordering
metadata <- metadata %>%
  dplyr::mutate(
    Diagnosis = factor(Diagnosis, levels = c("FTA", "minimallyFTC", "widelyFTC")),
    Tumor_type = factor(Tumor_type, levels = c("FTA", "FTC")),
    Metastasis = factor(Metastasis, levels = c("F", "T")),
    Ras = factor(Ras, levels = c("F", "T")),
    invasive = factor(invasive, levels = c("F", "T")),
    label = gsub("-", "_", ID),
    mtDNA_copynumber = round(mtDNA_copynumber),
    Tumor_size_factor = ifelse(Tumor_size < 1, "< 1cm",
                               ifelse(Tumor_size < 2, "< 2cm", "< 4cm")) %>% factor(.),
    Driver_mut_thyroid = Driver_mut,
    Patient = factor(Patient, levels = unique(patient_order))
  ) %>%
  dplyr::arrange(Diagnosis, Ras, Patient, Metastasis, invasive, Tip_total_mut, distance_from_invasive_lenient)

# FIGURE 2 ----------------------------------------------------------------
## extract information of median for each category (FTA, invasive of minimally FTC, noninvasive of minimally FTC, invasive of widely FTC, noninvasive of widely FTC) 
paired_metadata <- metadata %>% dplyr::filter(Tumor == "T" & merged_bam == "F" & High_depth == "T" & Metastasis == "F" ) %>% 
  dplyr::select(ID, Diagnosis, Patient, invasive, Tip_snv_mut, Diagnosis_age, Ras, cna_tip, Tip_total_mut, MRCA_snv_mut,
                Driver_mut_num, Driver_mut_num_cosmic, Tip_driver_mut_cosmic) %>%
  dplyr::mutate(
    Tip_total_driver_num = Tip_driver_mut_cosmic + Driver_mut_num,
    total_driver_num = Driver_mut_num + Driver_mut_num_cosmic, 
    tip_snv_mutation_rate = Tip_snv_mut / Diagnosis_age,
                tip_total_mutation_rate = Tip_total_mut / Diagnosis_age,
                mrca_snv_mutation_rate = MRCA_snv_mut / Diagnosis_age,
                branch_snv_mut = Tip_snv_mut - MRCA_snv_mut,
                branch_ratio = Tip_snv_mut / MRCA_snv_mut) 

paired_metadata$Patient <- factor(paired_metadata$Patient, levels = unique(paired_metadata$Patient))

avg_cn_dat <- lapply(paired_metadata$ID, function(ID){
  clean_segment_file <- list.files(path = file.path(DATA_DIR, "50_additional/clean_segments"), 
                                   pattern = paste0(ID,"_N.",".*clean$"), full.names = T)
  clean_seg <- read.table(clean_segment_file, header = T, comment.char = "", col.names = c('CHROM','START_POS','END_POS','TOTCN','MAJCN','MINCN')) %>%
    dplyr::mutate(cn_width = TOTCN * (END_POS - START_POS), width = END_POS - START_POS)
  avg_totcn <- sum(clean_seg[["cn_width"]]) / sum(clean_seg[["width"]])
  return(tibble(avg_cn = avg_totcn, ID = ID))
}) %>% data.table::rbindlist()

paired_metadata_cn <- left_join(x = paired_metadata, y = avg_cn_dat, by = 'ID')


# Dummy rows to add  ------------------------------------------------------

paired_metadata_cn <- paired_metadata_cn %>% dplyr::mutate(Diagnosis_spaced = Diagnosis)
# add a dummy row
dummy_row <- tibble::tibble(
  ID = NA, Diagnosis = NA, Patient = NA, invasive = NA, Tip_snv_mut = NA, Diagnosis_age = NA,
  Ras = NA, cna_tip = NA, Tip_total_mut = NA, MRCA_snv_mut = NA,
  Driver_mut_num = NA, Tip_driver_mut_cosmic = NA,
  Tip_total_driver_num = NA, tip_snv_mutation_rate = NA,
  tip_total_mutation_rate = NA, mrca_snv_mutation_rate = NA,
  branch_snv_mut = NA, branch_ratio = NA, avg_cn = NA,
  Diagnosis_spaced = factor("spacer1", levels = c("FTA", "spacer1", "minimallyFTC", "widelyFTC"))
)

# combine
paired_metadata_cn <- dplyr::bind_rows(paired_metadata_cn, dummy_row)

# x positions of the groups (used below and in the plots)
X_BREAKS = list(
  "FTA_vs_miFTC_vs_wiFTC" =  c(0.8, 1.6, 2.2),
  "FTC-body_vs_invasive" = c(1.1, 1.3, 1.5, 1.7)
)

# SUBFIGURE A. FTA VERSUS FTC ---------------------------------------------
patient_paired_metadata_cn_list <- list()

patient_paired_metadata_cn_list[['FTA_vs_miFTC_vs_wiFTC']] <- paired_metadata_cn %>% group_by(Patient, Diagnosis) %>%   # group by two keys
  mutate(Tip_snv_mut_mean = mean(Tip_snv_mut, na.rm = TRUE), 
         total_driver_num = total_driver_num, 
         #Tip_total_driver_num_mean = mean(Tip_total_driver_num, na.rm = TRUE),
         CN_mean = mean(avg_cn, na.rm = TRUE),
         Diagnosis = factor(Diagnosis, levels = c("FTA","minimallyFTC", "widelyFTC")),
         Diagnosis_spaced = factor(Diagnosis_spaced, levels = c('FTA','spacer1','minimallyFTC','widelyFTC')),
         xPos = case_when(
           Diagnosis == "FTA"           ~ X_BREAKS[[1]][[1]],   # position of the first group
           Diagnosis == "minimallyFTC"         ~ X_BREAKS[[1]][[2]], # position of the second group (reduce to 1.2 to place it closer to FTA)
           Diagnosis == "widelyFTC"         ~ X_BREAKS[[1]][[3]]    # position of the third group
         ))  %>%
  ungroup() %>% dplyr::distinct(Patient, Diagnosis, .keep_all = T)

patient_paired_metadata_cn_list[['FTC-body_vs_invasive']] <- paired_metadata_cn %>% dplyr::filter(Diagnosis == 'minimallyFTC' | Diagnosis == 'widelyFTC') %>% group_by(Diagnosis, Patient, invasive) %>%   # group by two keys
  mutate(Tip_snv_mut_mean = mean(Tip_snv_mut, na.rm = TRUE), 
         total_driver_num = total_driver_num, 
         Tip_total_driver_num_mean = mean(Tip_total_driver_num, na.rm = TRUE),
         CN_mean = mean(avg_cn, na.rm = TRUE),
         invasive = factor(invasive, levels = c("F","T")), 
         Diagnosis_invasive = paste0(Diagnosis,'_',invasive), 
         xPos = case_when(
           Diagnosis_invasive == 'minimallyFTC_F' ~ X_BREAKS[[2]][[1]],
           Diagnosis_invasive == 'minimallyFTC_T' ~ X_BREAKS[[2]][[2]],
           Diagnosis_invasive == 'widelyFTC_F' ~ X_BREAKS[[2]][[3]], 
           Diagnosis_invasive == 'widelyFTC_T' ~ X_BREAKS[[2]][[4]]
         ))  %>%
  ungroup() %>% dplyr::distinct(Patient, invasive, .keep_all = T)



# COMPARING P VALUES  -----------------------------------------------------

### Wilcoxon tests between diagnostic groups
# input data frame (e.g. patient_paired_metadata_cn_list[['FTA_vs_miFTC_vs_wiFTC']])
df <- patient_paired_metadata_cn_list[['FTA_vs_miFTC_vs_wiFTC']] %>% 
  mutate(Big_diagnosis = ifelse(grepl('FTC', Diagnosis), 'FTC', 'FTA'))

# drop rows with NA Diagnosis
df <- subset(df, !is.na(Diagnosis))

# Diagnosis as a factor in the order "FTA", "minimallyFTC", "widelyFTC"
df$Diagnosis <- factor(df$Diagnosis, levels = c("FTA", "minimallyFTC", "widelyFTC"))

# variables to compare
vars_to_test <- c("Tip_snv_mut", "total_driver_num", "avg_cn")

# Exact Wilcoxon rank-sum test (see header): uses the exact conditional distribution even with ties; base wilcox.test(exact = TRUE) falls back to the normal approximation when ties are present.
wx_exact <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  dd <- data.frame(y = c(a, b), g = factor(rep(c("a", "b"), c(length(a), length(b)))))
  as.numeric(coin::pvalue(coin::wilcox_test(y ~ g, data = dd, distribution = "exact")))
}
pairwise_exact <- function(x, g) {   # lower-triangular matrix in the same shape as pairwise.wilcox.test(p.adjust.method = "none")
  lv <- levels(droplevels(factor(g)))
  m <- matrix(NA_real_, length(lv) - 1, length(lv) - 1, dimnames = list(lv[-1], lv[-length(lv)]))
  for (i in 2:length(lv)) for (j in 1:(i - 1)) m[lv[i], lv[j]] <- wx_exact(x[g == lv[j]], x[g == lv[i]])
  m
}

# pairwise exact Wilcoxon test for each variable (no p-value adjustment)
for (var in vars_to_test) {
  cat("---------------------------------------------------\n")
  cat("Variable:", var, "\n")
  
  # pairwise Wilcoxon tests (no p-value adjustment)
  result <- pairwise_exact(df[[var]], df$Diagnosis)   # exact version of pairwise.wilcox.test(..., exact = FALSE)
  
  print(result)
  cat("\n")
  
  result2 <- pairwise_exact(df[[var]], df$Big_diagnosis)   # exact test (see header)
  print(result2)
  cat('\n')
}
## pvalue tip_total_driver_num
### FTA-miFTC : 0.67, FTA-wiFTC : 0.70, miFTC-wiFTC : 0.67   (exact test)
mini_total_driver_mean <- df %>% dplyr::filter(Diagnosis == 'minimallyFTC') %>% dplyr::pull(total_driver_num) %>% mean()
wide_total_driver_mean <- df %>% dplyr::filter(Diagnosis == 'widelyFTC') %>% dplyr::pull(total_driver_num) %>% mean()

FTA_total_driver_mean <- df %>% dplyr::filter(Big_diagnosis == 'FTA') %>% dplyr::pull(total_driver_num) %>% mean()
FTC_total_driver_mean <- df %>% dplyr::filter(Big_diagnosis == 'FTC') %>% dplyr::pull(total_driver_num) %>% mean()

#
## pvalue avg_cn
### FTA-miFTC:0.28, FTA-wiFTC: 0.20, miFTC-wiFTC: 0.044   (exact test)

mini_avg_cn <- df %>% dplyr::filter(Diagnosis == 'minimallyFTC') %>% dplyr::pull(avg_cn)
#[1] 1.916622 1.930626 2.001843 2.174540 2.834482 2.001393 2.044458 2.031906
wide_avg_cn <- df %>% dplyr::filter(Diagnosis == 'widelyFTC') %>% dplyr::pull(avg_cn)
# 3.216808 3.462386
fta_avg_cn <- df %>% dplyr::filter(Diagnosis == 'FTA') %>% dplyr::pull(avg_cn)
# 1.978141 2.008445 1.917660
wx_exact(mini_avg_cn, wide_avg_cn)   # exact version of wilcox.test(x, y)
wx_exact(fta_avg_cn, wide_avg_cn)
wx_exact(fta_avg_cn, mini_avg_cn)

### Wilcoxon tests: invasive vs. non-invasive
# variables to analyse
vars_to_test <- c("Tip_snv_mut", "Tip_total_driver_num", "avg_cn")

# data
patient_paired_metadata_cn_list[['FTA_vs_miFTC_vs_wiFTC']]
df <- patient_paired_metadata_cn_list[['FTC-body_vs_invasive']] %>%
  filter(Diagnosis %in% c("minimallyFTC", "widelyFTC"))

# function: wilcox.test by Diagnosis for one variable
get_group_pvals <- function(var) {
  df %>%
    group_by(Diagnosis) %>%
    summarise(
      variable = var,
      p_value = wx_exact(cur_data()[[var]][cur_data()$invasive == "F"], cur_data()[[var]][cur_data()$invasive == "T"]),   # exact test (see header)
      .groups = "drop"
    )
}

# run for each variable and combine the results
result <- map_dfr(vars_to_test, get_group_pvals)

result

# Running plots  ----------------------------------------------------------
## pvalue tip_total_driver_num
### FTA-miFTC : 0.67, FTA-wiFTC : 0.70, miFTC-wiFTC : 0.67   (exact test)
## pvalue avg_cn
### FTA-miFTC:0.28, FTA-wiFTC: 0.20, miFTC-wiFTC: 0.044   (exact test)

ANNOTATIONS <- 
  list(
    "FTA_vs_miFTC_vs_wiFTC" = list(
      # 3 comparisons => length 3 (same order as COMPARISONS[["FTA_vs_miFTC_vs_wiFTC"]])
      "DRIVER_mut" = c('0.67','0.70','0.67'),   # exact test
      "CN"         = c('0.28','0.20','0.04')    # exact test
    ),
    "FTC-body_vs_invasive" = list(
      # 1 comparison => length 1
      "DRIVER_mut" = c('0.997','1'),       # exact test
      "CN"         = c('0.574','0.667')    # exact test
    )
  )

# check that the p-values printed on the figure match the exact-test values above; stop otherwise (both panels)
local({
  d <- patient_paired_metadata_cn_list[['FTA_vs_miFTC_vs_wiFTC']] %>% dplyr::filter(!is.na(Diagnosis))
  d$Diagnosis <- factor(d$Diagnosis, levels = c("FTA", "minimallyFTC", "widelyFTC"))
  pw <- function(v) {
    m <- pairwise_exact(d[[v]], d$Diagnosis)   # exact test (see header)
    sprintf('%.2f', c(m["minimallyFTC", "FTA"], m["widelyFTC", "FTA"], m["widelyFTC", "minimallyFTC"]))  # order of COMPARISONS
  }
  stopifnot(identical(pw("total_driver_num"), ANNOTATIONS[["FTA_vs_miFTC_vs_wiFTC"]][["DRIVER_mut"]]),
            identical(pw("avg_cn"),           ANNOTATIONS[["FTA_vs_miFTC_vs_wiFTC"]][["CN"]]))
  b <- patient_paired_metadata_cn_list[['FTC-body_vs_invasive']] %>% dplyr::filter(Diagnosis %in% c("minimallyFTC", "widelyFTC"))
  pb <- function(v) sapply(c("minimallyFTC", "widelyFTC"), function(dg) {
    x <- b[b$Diagnosis == dg, ]; as.character(signif(wx_exact(x[[v]][x$invasive == "F"], x[[v]][x$invasive == "T"]), 3))
  }, USE.NAMES = FALSE)
  stopifnot(identical(pb("Tip_total_driver_num"), ANNOTATIONS[["FTC-body_vs_invasive"]][["DRIVER_mut"]]),
            identical(pb("avg_cn"),               ANNOTATIONS[["FTC-body_vs_invasive"]][["CN"]]))
})

color_list <- list(
  "FTA_vs_miFTC_vs_wiFTC" = c("FTA" = patient_cols[1], "minimallyFTC" = patient_cols[2] , "widelyFTC" = patient_cols[3]),
  "FTC-body_vs_invasive" = c("minimallyFTC_F" = patient_cols[2], "minimallyFTC_T" = patient_cols[2] ,
                             "widelyFTC_F" = patient_cols[3], "widelyFTC_T" = patient_cols[3])
)
shape_list <- list(
  "FTA_vs_miFTC_vs_wiFTC" = c("FTA" = 16, "minimallyFTC" = 16, "widelyFTC" = 16),
  "FTC-body_vs_invasive" = c("F" = 16, "T" = 16)
)


X_VAR = list('FTA_vs_miFTC_vs_wiFTC' = "Diagnosis_spaced",
             'FTC-body_vs_invasive' = "Diagnosis"
)

SHAPE_VAR = list('FTA_vs_miFTC_vs_wiFTC' = "Diagnosis_spaced",
                 'FTC-body_vs_invasive' = "invasive"
)

X_LABELS <- list('FTA_vs_miFTC_vs_wiFTC' = c('FTA','','miFTC','wiFTC'),
                 'FTC-body_vs_invasive' = c('miFTC','wiFTC')
)
X_MEZZLABEL <- c('minimally FTC', 'widely FTC')
X_SUPERLABEL <- c('FTC')
COMPARISONS <- list(
  'FTA_vs_miFTC_vs_wiFTC' = list(c('FTA','minimallyFTC'), c('FTA','widelyFTC'), c('minimallyFTC','widelyFTC')),
  'FTC-body_vs_invasive' = list(c('F','T'))
)

AXIS_TEXT_SIZE = 8
AXIS_TITLE_SIZE = 8
SIGNIF_TEXT_SIZE = 2

BOXPLOT_WIDTH = 0.4
JITTER_WIDTH = 0.1
JITTER_SIZE = 0.6

SNV_mut_MAX = 3200
DRIVER_mut_MAX = 5
CN_MAX = 5

X_POSITIONS_MIN <- list(
  "FTA_vs_miFTC_vs_wiFTC" =  c(X_BREAKS[[1]][1],X_BREAKS[[1]][1],X_BREAKS[[1]][2]),
  "FTC-body_vs_invasive" = c(X_BREAKS[[2]][1],X_BREAKS[[2]][3])
)
X_POSITIONS_MAX <- list(
  "FTA_vs_miFTC_vs_wiFTC" =  c(X_BREAKS[[1]][2],X_BREAKS[[1]][3],X_BREAKS[[1]][3]),
  "FTC-body_vs_invasive" = c(X_BREAKS[[2]][2],X_BREAKS[[2]][4])
)
X_LABELS = list(
  "FTA_vs_miFTC_vs_wiFTC" = c("FTA", "miFTC", "wiFTC"),
  "FTC-body_vs_invasive" = c('Non-invasive','Invasive','Non-invasive','Invasive')
)
FILL_VAR = list(
  "FTA_vs_miFTC_vs_wiFTC" = 'Diagnosis_spaced',
  "FTC-body_vs_invasive" = 'Diagnosis_invasive'
)
Y_POSITIONS <- list(
  "FTA_vs_miFTC_vs_wiFTC" = list(
    # 3 comparisons => length 3 (same order as COMPARISONS[["FTA_vs_miFTC_vs_wiFTC"]])
  #  "SNV_mut"    = c(SNV_mut_MAX - 800, SNV_mut_MAX - 500, SNV_mut_MAX - 200),
    "DRIVER_mut" = c(DRIVER_mut_MAX - 1, DRIVER_mut_MAX - 0.6, DRIVER_mut_MAX - 0.2),
    "CN"         = c(CN_MAX - 1, CN_MAX - 0.6, CN_MAX - 0.2)
  ),
  "FTC-body_vs_invasive" = list(
    # 1 comparison => length 1
   # "SNV_mut"    = c(SNV_mut_MAX - 200,SNV_mut_MAX - 200),  # 1 comparison => length 1
    "DRIVER_mut" = c(DRIVER_mut_MAX - 0.2, DRIVER_mut_MAX - 0.2),
    "CN"         = c(CN_MAX -0.2, CN_MAX -0.2)
  )
)

Y_LIMITS = list(#'SNV_mut' = c(0,SNV_mut_MAX), 
                'DRIVER_mut' = c(0,DRIVER_mut_MAX), 
                'CN' = c(0,CN_MAX))
Y_LABELS = list(#'SNV_mut' = c(0,500,1000,1500,2000), 
                'DRIVER_mut' = c(0,1,2,3,4), 
                'CN' = c(0,1,2,3,4))
Y_LAB = list(#'SNV_mut' = 'Number of base substitutions',
             'DRIVER_mut' = 'Number of point driver mutations',
             'CN' = 'Average copy number')


# RUNNING FOR THE FOR LOOP
# comparison groups (rows)
row_keys <- c("FTA_vs_miFTC_vs_wiFTC", 
              "FTC-body_vs_invasive")

# metrics to plot
metrics <- c( "DRIVER_mut", "CN") # "SNV_mut",

# variable name in the data set for each metric
y_vars <- list(
  #SNV_mut = "Tip_snv_mut_mean",
  DRIVER_mut = "total_driver_num", # Tip_total_driver_num_mean
  CN = "CN_mean"
)


all_plots <- list()
index <- 1
for (rk in row_keys) {
  for (m in metrics) {
    y_var <- y_vars[[m]]
    if (m == 'DRIVER_mut' & rk == 'FTC-body_vs_invasive'){
      y_var <- "Tip_total_driver_num_mean"
    } 
    
    signif_df <- tibble(
      xmin = X_POSITIONS_MIN[[rk]],
      xmax = X_POSITIONS_MAX[[rk]],
      annotations = ANNOTATIONS[[rk]][[m]],
      y_position = Y_POSITIONS[[rk]][[m]]
    ) %>% dplyr::mutate(group = row_number())
    
    if (rk =='FTA_vs_miFTC_vs_wiFTC'){
      y2lns <- -0.8
      y2lbls <- y2lns - 0.2
      x2lns <- data.frame(xstart = c(X_BREAKS[[2]][1], X_BREAKS[[2]][3]), xend = c(X_BREAKS[[2]][2],X_BREAKS[[2]][4]))
      x2lbls <- data.frame(x = c( (X_BREAKS[[2]][1] + X_BREAKS[[2]][2]) / 2, (X_BREAKS[[2]][3] + X_BREAKS[[2]][4])/2 ), lbl = c('minimally', 'widely'))
      y1lns <- y2lns - 0.5
      y1lbls <- y1lns -0.2
      x1lns <- data.frame(xstart = X_BREAKS[[1]][2] , xend = X_BREAKS[[1]][3] )
      x1lbls <- data.frame(x= (X_BREAKS[[1]][2] + X_BREAKS[[1]][3] )/ 2, lbl = c('FTC'))
      LINETHICKENSS = 0.25
      
      summarised_data <- patient_paired_metadata_cn_list[[rk]] %>%
        group_by(xPos) %>%
        summarise(mean_y = mean(!!sym(y_var)),
                  se_y = sd(!!sym(y_var)) / sqrt(n()), 
                  Diagnosis_spaced = unique(Diagnosis_spaced))
      
      p <- ggplot(summarised_data, aes(x = xPos, y = mean_y)) +
        geom_bar(stat = "identity", aes_string(fill = FILL_VAR[[rk]])) +
        geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), width = 0.2) +
        geom_jitter(data = patient_paired_metadata_cn_list[[rk]], 
          aes_string(x = "xPos", y = y_var, fill = FILL_VAR[[rk]]), colour= 'black', pch = 21, #, shape = SHAPE_VAR[[rk]] #X_VAR[[rk]]
          position = position_jitter(width = JITTER_WIDTH),
          size = JITTER_SIZE,
          alpha = 1
        ) +       # font / axis settings
        coord_cartesian(ylim=c(-0.2,NA), clip="off") +
        theme_classic()  + 
        theme(
          axis.text = element_text(size = AXIS_TEXT_SIZE),
          axis.title = element_text(size = AXIS_TITLE_SIZE),
          axis.title.x=element_blank(),
          plot.margin=margin(b=0.5,t=0.5,l=0.6,r=0.6,unit="lines"),
          base_family = 'Arial'
        ) + 
        
        scale_x_continuous(
          breaks = X_BREAKS[[rk]], # positions shown on the x axis
          labels = X_LABELS[[rk]]  # labels for those positions
        ) +
        
        geom_segment(data=x1lns,mapping=aes(x=xstart,xend=xend), y = y1lns, yend = y1lns, #, y= ystart, yend = yend), 
                     linewidth=LINETHICKENSS, inherit.aes = FALSE) +
        geom_text(data=x1lbls,mapping=aes(x=x,label=lbl, y = y), y = y1lbls, 
                  size=8/.pt,  inherit.aes = FALSE) + 
        
        # manual colour / shape scales
        scale_color_manual(values = color_list[[rk]]) +
        scale_fill_manual(values = color_list[[rk]]) +
        scale_shape_manual(values = shape_list[[rk]]) +
        
        # y-axis range / labels
        scale_y_continuous(
          expand = expansion(mult = c(0.05, 0.05)),
          limits = Y_LIMITS[[m]],
          breaks = Y_LABELS[[m]],
          labels = Y_LABELS[[m]]
        ) +
        ylab(Y_LAB[[m]]) +
        xlab('') +
        
        # comparison brackets
        ggpubr::geom_signif(
          data = signif_df,
          aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position, group = group),
          manual = TRUE,
          textsize = SIGNIF_TEXT_SIZE,
          family = 'Arial',
          inherit.aes = FALSE, 
          color = 'black'
        ) + theme(legend.position="none")
      
    } else {
      y2lns <- -0.8
      y2lbls <- y2lns - 0.2
      x2lns <- data.frame(xstart = c(X_BREAKS[[2]][1], X_BREAKS[[2]][3]), xend = c(X_BREAKS[[2]][2],X_BREAKS[[2]][4]))
      x2lbls <- data.frame(x = c( (X_BREAKS[[2]][1] + X_BREAKS[[2]][2]) / 2, (X_BREAKS[[2]][3] + X_BREAKS[[2]][4])/2 ), lbl = c('minimally', 'widely'))
      #y1lns <- y2lns - 0.5
      #y1lbls <- y1lns -0.2
      #x1lns <- data.frame(xstart = (X_BREAKS[[2]][1] + X_BREAKS[[2]][2]) / 2, xend = (X_BREAKS[[2]][3] + X_BREAKS[[2]][4])/2 )
      #x1lbls <- data.frame(x= (X_BREAKS[[2]][1] + X_BREAKS[[2]][2] + X_BREAKS[[2]][3] + X_BREAKS[[2]][4])/ 4, lbl = c('FTC'))
      LINETHICKENSS = 0.25
      p <- ggplot(patient_paired_metadata_cn_list[[rk]], 
                  aes_string(x = "xPos", y = y_var,  # X_VAR[[rk]]
                             color = FILL_VAR[[rk]])) +  #
        geom_line(aes_string(x = "xPos", y = y_var, group = "Patient", colour = FILL_VAR[[rk]])) + # X_VAR[[rk]]
        geom_jitter(
          aes_string(x = "xPos", y = y_var, fill = FILL_VAR[[rk]]), colour= 'black', pch = 21, #, shape = SHAPE_VAR[[rk]] #X_VAR[[rk]]
          position = position_jitter(width = 0, height = 0),
          size = JITTER_SIZE,
          alpha = 1
        ) + 
        coord_cartesian(ylim=c(0,NA), clip="off") +
        theme_classic()  + 
        theme(axis.title = element_text(size = AXIS_TITLE_SIZE),
              axis.title.x=element_blank(),
              axis.text = element_text(size = AXIS_TEXT_SIZE),
              axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = 0, b = 10)),  # 10 pt margin below, 0 pt above
              plot.margin=margin(b=0.5,t=0.5,l=0.6,r=0.6,unit="lines"),
              base_family = 'Arial'
        ) +
        geom_segment(data=x2lns,mapping=aes(x=xstart,xend=xend), y = y2lns, yend = y2lns, #, y= ystart, yend = yend), 
                     linewidth=LINETHICKENSS, inherit.aes = FALSE) +
        geom_text(data=x2lbls,mapping=aes(x=x,label=lbl, y = y), y = y2lbls, 
                  size=8/.pt, inherit.aes = FALSE) + 
        scale_x_continuous(
          limits = c(min(X_BREAKS[[rk]]) - 0.05, max(X_BREAKS[[rk]]) + 0.05),
          breaks = X_BREAKS[[rk]],               # positions shown on the x axis
          labels = X_LABELS[[rk]]  # labels for those positions
        ) +
        # manual colour / shape scales
        scale_color_manual(values = color_list[[rk]]) +
        scale_fill_manual(values = color_list[[rk]]) +
        scale_shape_manual(values = shape_list[[rk]]) +
        
        # y-axis range / labels
        scale_y_continuous(
          limits = Y_LIMITS[[m]],
          breaks = Y_LABELS[[m]],
          labels = Y_LABELS[[m]]
        ) +
        ylab(Y_LAB[[m]]) +
        xlab('') +
        
        # comparison brackets
        ggpubr::geom_signif(
          data = signif_df,
          aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position, group = group),
          manual = TRUE,
          textsize = SIGNIF_TEXT_SIZE,
          family = 'Arial',
          inherit.aes = FALSE, 
          color = 'black'
        ) + theme(legend.position="none")
    }
    all_plots[[index]] <- p #+ theme(plot.margin = margin(t=10, r=5, b=10, l=5))
    index <- index + 1
  }
}

pdf('Fig2b-e_drivers_and_ploidy.pdf', width = 7.5, height = 3)
ggarrange(align = 'h', plotlist = all_plots, nrow = 1)
dev.off()
