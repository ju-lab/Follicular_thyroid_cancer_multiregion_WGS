# author: taewoo
# date: Oct 18 2025 
# description: Fig. 4a (clock-like mutation rate), Fig. 4c (timing of copy-number gains) and
#              Fig. 4d (tumor expansion time vs. tumor size)

# Paths: run from the repository root; input data are read from data/ and figures are written to figures/.
# Fig. 4c uses combine_axes() from github.com/davidearn/plague_growth (GPL-3.0), sourced at a fixed
# commit, so this script needs internet access.
CODE_DIR   <- normalizePath(".")
DATA_DIR   <- file.path(CODE_DIR, "data")
FIGUREPATH <- file.path(CODE_DIR, "figures")
dir.create(FIGUREPATH, showWarnings = FALSE)
source(file.path(CODE_DIR, 'Thyroid_sources.R'))
setwd(FIGUREPATH)

# (FIGURE 4) --------------------------------------------------------------
# (FIGURE 4.A) Mutation per Year with samples more than 5  -----------------------------
metadata_file = file.path(DATA_DIR, "thyroid_sample_metadata.xlsx")
metadata <- read_sample_metadata(metadata_file)

# patient order (by Diagnosis + Ras)
patient_order <- metadata %>% 
  dplyr::filter(invasive == 'F') %>% 
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
    specific_type = paste0(Diagnosis, "_", invasive),
    specific_ras = paste0(Ras, "-", Diagnosis),
    Patient = factor(Patient, levels = unique(patient_order))
  ) %>%
  dplyr::arrange(Diagnosis, Ras, Patient, Metastasis, invasive, Tip_total_mut, distance_from_invasive_lenient)


mutation_dat <-  metadata %>% 
  dplyr::filter(Tumor == "T" & merged_bam == "F" & High_depth == "T") %>%
  dplyr::group_by(Patient) %>% 
  dplyr::reframe(n = n(), MRCA_snv_mut = MRCA_snv_mut, Diagnosis_age = Diagnosis_age, Diagnosis = Diagnosis) %>% 
  dplyr::distinct(Patient, .keep_all = T) %>%
  dplyr::mutate(cutoff = n>=5)

# calculate the regression line 
# 50 per year for MRCA mutations, 1904 mutations ~ 38 years
mutation.lm.all <- lm(data = mutation_dat %>% dplyr::filter(n >= 5), formula = MRCA_snv_mut ~ Diagnosis_age)


# correlation after filtering n >= 5
df <- mutation_dat %>% filter(n >= 5)

ct <- cor.test(df$Diagnosis_age, df$MRCA_snv_mut, method = "pearson")
r <- unname(ct$estimate)
p <- ct$p.value

# print the result
cat(sprintf("r = %.3f, p = %.4f\n", r, p))

# calculate the confidence interval of regression line slope
## first method
#https://stats.stackexchange.com/questions/221630/how-to-calculate-the-confidence-interval-of-the-x-intercept-in-a-linear-regressi
dat <- mutation_dat %>% dplyr::filter(n>=5) %>% dplyr::select(MRCA_snv_mut, Diagnosis_age)%>%
  dplyr::rename(y = MRCA_snv_mut, x = Diagnosis_age)

xInterceptCI <- function(x, alpha = 0.05, ...) {
  intercept <- coef(x)[1]              
  slope <- coef(x)[2]                  
  meanX <- mean(x$model[,2])           
  n <- length(x$model[,2])             
  tstar <- qt(alpha/2, n-2)            
  sxx <- sum(x$model[,2]^2) - sum(x$model[,2])^2 / n
  SSresidual <- (1-cor(x$model[,1], x$model[,2])^2) * 
    (sum(x$model[,1]^2)-sum(x$model[,1])^2/n)
  S <- sqrt(SSresidual/(n-2))
  SEslope <- S / sqrt(sxx)
  Xintercept <- - intercept / slope
  y0 <- 0
  g <- (tstar / (slope/SEslope))^2
  left <- (Xintercept - meanX) * g
  bottom <- 1 - g
  Right <- (tstar * S / slope) * sqrt( ((Xintercept - meanX)^2/sxx) + bottom/n)
  lower <- Xintercept + (left + Right) / bottom
  upper <- Xintercept + (left - Right) / bottom
  return(c(lower,upper))
}

xInterceptCI(mutation.lm.all,alpha = 0.05) # 15.05835, 45.68675

## decorate lm object with a new class lm_left (the fitted line is drawn only up to the largest x)
lm_left <- function(formula,data,...){
  mod <- lm(formula,data)
  class(mod) <- c('lm_left',class(mod))
  mod
}

predictdf.lm_left <-
  function(model, xseq, se, level){
    init_range = range(model$model$x)
    ## here the main code: truncate to x values at the left
    xseq <- xseq[xseq <=init_range[2]]
    ggplot2:::predictdf.default(model, xseq[-length(xseq)], se, level)
  }

##I referred
##https://stackoverflow.com/questions/26705554/extend-geom-smooth-in-a-single-direction
mutation_age_plot <- ggplot(mutation_dat, aes(x = Diagnosis_age, y = MRCA_snv_mut)) +
  # confidence band
  stat_smooth(
    data = mutation_dat %>% dplyr::filter(n >= 5),
    method = "lm_left",
    geom = "ribbon",
    fullrange = TRUE,
    fill = "gray80",
    alpha = 0.5,
    se = TRUE
  ) +
  
  # line
  stat_smooth(
    data = mutation_dat %>% dplyr::filter(n >= 5),
    method = "lm_left",
    geom = "line",
    fullrange = TRUE,
    col = 'black',
    size = 0.5,
    se = FALSE
  ) +
  
  # points
  geom_point(
    data = mutation_dat,
    aes(colour = Diagnosis, size = ifelse(n>=5, n, 1)) #, alpha = cutoff
  ) +
  scale_color_manual(
    values = c("FTA" = patient_cols[1], "minimallyFTC" = patient_cols[2], "widelyFTC" = patient_cols[3]),
    labels = c('FTA','Minimally FTC','Widely FTC')
  ) +
  scale_alpha_manual(
    name = "Patient with more than 5 samples",
    values = c("TRUE" = 1, "FALSE" = 0.4),
    labels = c('TRUE', 'FALSE')
  ) +
  scale_size_continuous(
    name = 'No. samples',
    breaks = c(1, 5, 9),
    labels = c('< 5', 5, 9),
    range = c(1, 5)  # smaller maximum point size
  ) +
  scale_x_continuous(limits = c(0, 80)) +
  coord_cartesian(ylim = c(0, 2500)) +
  theme_classic() +
  ylab("Number of base substitutions in MRCA") +
  xlab("Age at diagnosis") +
  guides(color = guide_legend(order = 1, override.aes = list(size = 2)),
         size  = guide_legend(order = 2),
         alpha = 'none'
         ) + 
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7)
  ) + 
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2,0.7),
    legend.frame = NULL, 
    legend.direction = "vertical",
    legend.box = "vertical"
  )

# (FIGURE 4.C) timing of copy-number gains ------------------------------------
source(file.path(CODE_DIR, 'TIMINGANALYSIS_AMPONLY.R'))
setwd(FIGUREPATH)  # TIMINGANALYSIS_AMPONLY.R changes the working directory
final_annot

copy_timing <- read_excel(metadata_file, sheet = "copynumber_timing2") %>%
  filter(Patient %in% c('C08','C10','C04','C07','C01','C09','C02','C03'))  %>% 
  dplyr::mutate(year = Mutations / mutation.lm.all$coefficients[["Diagnosis_age"]],
                tip_year = Mutations / mutation.lm.all$coefficients[["Diagnosis_age"]])  %>% 
  rowwise() %>% 
  dplyr::mutate(RAS = ifelse(Patient %in% c("C07","C08","C06","C10"), "T","F" ),
                RAS = factor(RAS, levels = c('F','T')),
                Diagnosis = ifelse(Patient %in% c('C10','C09'), 'widelyFTC', 'minimallyFTC'),
                Patient = factor(Patient, levels = c('C08','C10','C04','C07','C01','C09','C02','C03') |> rev()),
                cna_type = ifelse(grepl('WCA',cna_type), 'WGD', ifelse(is.na(cna_type), 'No event',cna_type)),
                cna_type = cna_type %>% gsub(';',' & ',.) %>% gsub('_', ' ',. ),
                Event = factor(Event, levels = c('MRCA','branch'))) %>%
  dplyr::arrange(RAS, Patient, year) %>% 
  dplyr::mutate(cna_type = factor(cna_type))

### ADDITION OF MRCA AND BRANCH EVENTS OF THE AMPLIFICATION 
### ADDING AMPLIFIED CHROMOSOMES AND THE TIMING IN THE PLOT 

coef_age <- as.numeric(mutation.lm.all$coefficients[['Diagnosis_age']])

# x-axis order of the timing panel (patients by age at diagnosis)
patient_levels <- copy_timing %>%
  filter((Patient == 'C10' | cutoff == "T") & Event == 'MRCA') %>%
  arrange(desc(Diagnosis_age)) %>% pull(Patient) %>% unique()

gap_df <- copy_timing %>%
  filter((Patient == 'C10' | cutoff == "T") & Event == 'MRCA') %>%
  mutate(
    Patient = factor(Patient, levels = patient_levels),
    Driver  = factor(ifelse(Driver == "RAS", "RAS", "Non-RAS"),
                     levels = c("RAS","Non-RAS")),
    y_pred  = Diagnosis_age * coef_age,
    y_low   = pmin(Mutations, y_pred),
    y_high  = pmax(Mutations, y_pred),
    y_mid   = (y_low + y_high) / 2,                 # tile centre
    height  = pmax(y_high - y_low, 0),              # tile height
    gap_type = if_else(Mutations >= y_pred, "above", "below")
  )

# add only the layer filling the gap between bar and prediction to g1 (placed after the bars)
g1 <- copy_timing %>%
  filter((Patient == 'C10' |  cutoff == "T") & Event == 'MRCA') %>%
  mutate(
    Patient = factor(Patient, levels = patient_levels),
    Driver  = factor(ifelse(Driver == "RAS","RAS","Non-RAS"),
                     levels = c("RAS","Non-RAS"))
  ) %>%
  ggplot() +
  geom_bar(aes(x = Patient, y = Mutations), fill = "black",
           width = 0.6, stat = "identity") +
  # fill the gap: aligned automatically on the discrete axis
  geom_tile(data = gap_df,
            aes(x = Patient, y = y_mid, height = height, width = 0.6,
                fill = gap_type),
            alpha = 0.28, color = NA, inherit.aes = FALSE) +
  geom_point(aes(x = Patient, y = Mutations), size = 1, color = "#D55E00") +
  geom_point(aes(x = Patient, y = Diagnosis_age * coef_age), size = 1) +
  scale_fill_manual(values = c(above = "#F0A35A", below = "#5AA0D6"),
                    name = "gap") +
  scale_y_continuous(sec.axis = dup_axis(name = NULL), limits = c(0, 4000)) +
  facet_grid(~ Driver, scales = "free", space = "free") +
  theme_classic() +
  theme(
    axis.text.y.left  = element_text(angle = 90, hjust = 0.4),
    axis.text.y.right = element_text(angle = 90, hjust = 0.4),
    axis.text.x       = element_text(angle=90),
    axis.title.y      = element_text(angle = 90, vjust = 0.5, size = 7),  # rotate y-axis title to vertical
    text              = element_text(family = "Arial"),
    axis.text         = element_text(size=7),
    axis.title        = element_text(size=7),
    plot.title        = element_text(size=7),
    legend.text       = element_text(size=7),
    legend.title      = element_text(size=7),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal"
  )


# 3) label data (same levels)
label_df <- final_annot %>%
  rename(Patient = patient) %>%
  mutate(Patient = factor(Patient, levels = patient_levels)) %>%
  left_join(
    copy_timing %>%
      filter((Patient == 'C10' | cutoff == "T") & Event == 'MRCA') %>%
      transmute(Patient = factor(Patient, levels = patient_levels),
                Driver  = factor(ifelse(Driver == "RAS", "RAS", "Non-RAS"),
                                 levels = c("RAS","Non-RAS"))),
    by = "Patient"
  ) %>%
  mutate(y_label = pmin(molecular_time_mean *
                          as.numeric(mutation.lm.all$coefficients[["Diagnosis_age"]]),
                        3950))

# 4) horizontal marker lines + vertical text (fixes alignment)
g1 <- g1 +
  geom_errorbar(
    data = label_df,
    aes(x = Patient, ymin = y_label, ymax = y_label, color = event),
    width = 0.55,              # length of the horizontal line (matched to bar width ~0.6)
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_df,
    aes(x = Patient, y = y_label +10, label = broad_amp, color = event),
    angle = 90, hjust = 0, vjust = 0.5, size = 7/.pt,
    inherit.aes = FALSE
  ) +
  scale_color_manual(name = 'Event', values = c(MRCA = "#D55E00", branch = "#0072B2"))



g2 <- copy_timing %>% 
  dplyr::filter( (Patient == 'C10' |  cutoff == "T") & Event == 'MRCA') |> 
  mutate(Patient = factor(Patient, levels = unique(copy_timing$Patient)), 
         Driver = case_when(
           Driver == "NA" ~ 'Non-RAS',
           Driver == 'RAS' ~ 'RAS'
         ), 
         Driver = factor(Driver, levels = c('RAS','Non-RAS'))) |> 
  ggplot() + geom_point(aes(x = Patient, Diagnosis_age), size = 1) + 
  theme_classic() + 
  labs(y = 'Age estimates based on\nclock-like mutation rate(50/yr)') + 
  scale_y_continuous(sec.axis = dup_axis(), limits = c(0,80)) +  #sec.axis = dup_axis(), 
  geom_bar(aes(x =Patient, year, fill = Diagnosis), position = 'dodge', stat = 'identity') + 
  theme(axis.text.y.left = element_text(angle = 90, hjust = 0.4),
        axis.text.y.right = element_text(angle = 90, hjust = 0.4),
        axis.title.y.right = element_text(angle= 90, hjust = 1),
        text = element_text(family = "Arial"),
        axis.text=element_text(size=7), #change font size of axis text
        axis.title=element_text(size=7), #change font size of axis titles
        plot.title=element_text(size=7), #change font size of plot title
        legend.text=element_text(size=7), #change font size of legend text
        legend.title=element_text(size=7),
        legend.position = NULL  # added
  )
#coord_flip()



# combine_axes(): David Earn's plague_growth repository (GPL-3.0), pinned to commit fc827bf
library(cowplot); library(gtable)
source("https://raw.githubusercontent.com/davidearn/plague_growth/fc827bfcb8edd812d508abb8b8662d6420998d31/analysis/plots/3axes.R")

copy_number_age_plot <-  combine_axes(g1,g2,add_pos="r", pad_inner = 0.5)

ggsave(plot = copy_number_age_plot, filename = 'Fig4c_copy_number_gain_timing.pdf', path = FIGUREPATH, width = 2.5, height = 3.5)


# (FIGURE 4.D) tumor expansion time vs. tumor size ------------------------------


copy_timing_subset <- copy_timing %>% 
  dplyr::filter( (Patient == 'C10' |  cutoff == "T") & Event == 'MRCA')
copy_timing_meta <- metadata |> distinct(Patient, .keep_all = T) |> filter(!is.na(Patient))
copy_timing_merge <- merge(copy_timing_subset, copy_timing_meta, by = 'Patient')





# tumor expansion time = age at diagnosis - age at the MRCA amplification
cor.test(copy_timing_merge$Diagnosis_age.x - copy_timing_merge$year, copy_timing_merge$Tumor_size, method = "pearson")

# linear fits: expansion time ~ tumor size, and tumor size ~ expansion time (growth rate)
tumor_expand_lm <- lm((Diagnosis_age.x - year) ~ Tumor_size, data = copy_timing_merge)
copy_timing_merge_flip  <- copy_timing_merge |> mutate(gap_time = Diagnosis_age.x - year)
tumor_expand_lm_flipped <- lm(Tumor_size  ~ (gap_time), data = copy_timing_merge_flip)
summary(tumor_expand_lm_flipped)

gap_timing <- ggplot(copy_timing_merge, aes(Tumor_size, Diagnosis_age.x - year)) +
    # confidence band
    stat_smooth(
      data = copy_timing_merge,
      method = "lm_left",
      geom = "ribbon",
      fullrange = TRUE,
      fill = "gray80",
      alpha = 0.5,
      se = TRUE
    ) +
    
    # line
    stat_smooth(
      data = copy_timing_merge,
      method = "lm_left",
      geom = "line",
      fullrange = TRUE,
      col = 'black',
      size = 0.5,
      se = FALSE
    ) +
  geom_point(size = 2, aes(color = Diagnosis.x, shape = Driver)) + 
  # dashed line from the origin: marks the lag phase before tumor expansion (drawn coordinates)
  geom_line(data = data.frame(x = c(0,1), y = c(0,32)), linetype = 2, aes(x,y)) + 
  theme_classic(base_family = 'Arial') +
  theme(axis.text=element_text(size=8), #change font size of axis text
        axis.title=element_text(size = 8), #change font size of axis titles
        plot.title=element_text(size=8), #change font size of plot title
        legend.text=element_text(size=7), #change font size of legend text
        legend.title=element_text(size=7),
        legend.margin = margin()
  ) +
  scale_color_manual(name = 'Diagnosis', values = c('minimallyFTC' =  patient_cols[2], 'widelyFTC' = patient_cols[3]),
                     labels = c('minimally FTC', 'widely FTC')) +
  scale_shape(name = 'RAS mutation', labels = c('RAS (-)', 'RAS (+)')) + 
  labs(x = 'Tumor size (cm)', 
       y = 'Tumor expansion time (yrs)') + 
  scale_y_continuous(limits = c(0, 50), labels = seq(0,50,10), breaks = seq(0,50,10)) + 
  scale_x_continuous(limits = c(-0.1,9.5), labels = c(0,2,4,6,8), breaks= c(0,2,4,6,8))

gap_timing_flipped <- gap_timing + coord_flip()










# prediction + 95% CI at Tumor_size = 0
pred0 <- predict(tumor_expand_lm, 
                 newdata = data.frame(Tumor_size = 0), 
                 interval = "confidence", 
                 level = 0.95)

pred0  # x-intercept of Fig. 4d (expansion time at tumor size 0) with 95% CI

# Save figures --------------------------------------------------------------
pdf('Fig4a_clock_mutation_rate.pdf', width = 3.5, height = 3.5)
mutation_age_plot
dev.off()

pdf('Fig4d_tumor_expansion_time.pdf', width = 3.1, height = 2.5)
gap_timing_flipped
dev.off()
