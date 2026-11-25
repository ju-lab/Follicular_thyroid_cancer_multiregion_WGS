# author: taewoo
# date: Oct 18 2025 
# description: Fig. 3a (phylogenetic trees), Fig. 3d (branch-to-trunk ratio by RAS status) and
#              Fig. 3e (mutation distance vs. physical distance between sister regions)

# Paths: run from the repository root; input data are read from data/ and figures are written to figures/.
CODE_DIR   <- normalizePath(".")
DATA_DIR   <- file.path(CODE_DIR, "data")
FIGUREPATH <- file.path(CODE_DIR, "figures")
dir.create(FIGUREPATH, showWarnings = FALSE)
source(file.path(CODE_DIR, 'Thyroid_sources.R'))
setwd(FIGUREPATH)

library(ggrepel)

# (FIGURE 3.a) Individual tree with HEATMAP  --------------------------------------------------------

tree_gg_list_final <- readRDS(file.path(DATA_DIR, "50_additional/sample_mutation_with_drivers/tree_gg_list.rds"))
## ADD A01 and A02 node color (#MAR 22 2025 TO ADD THE FTA TREE)
tree_gg_list_final_2 <- tree_gg_list_final
tree_gg_list_final_2$A01$tree <- tree_gg_list_final_2$A01$tree %>% dplyr::mutate(
  node_color = ifelse(label %in% c('A01_7_8','A01_3_4'), 1,
                      ifelse(label %in% c('A01_1_2','A01_5_6'), 2, NA)
  ) %>% factor(., levels = c(1,2,3,4)),
  label2 = c('A011','A012','A013','A014',NA,NA,NA) 
)
tree_gg_list_final_2$A03$tree <- tree_gg_list_final_2$A03$tree %>% dplyr::mutate(
  label2 = c('A031',NA,NA),
  node_color = c('9',NA,NA) %>% factor()
)
tree_gg_list_final_2$A02$tree <- tree_gg_list_final_2$A02$tree %>% dplyr::mutate(
  node_color = ifelse(label %in% c('A02_5_6','A02_3_4'), 1,
                      ifelse(label %in% c('A02_7_8','A02'), 2, NA)
  ) %>% factor(., levels = c(1,2,3,4)),
  label2 = c('A024','A021','A023','A022', NA, NA,NA)
)
tree_gg_list_final_2$C02$tree <- tree_gg_list_final_2$C02$tree %>% dplyr::mutate(
  driver = ifelse(!is.na(driver), 'TSHR;NFE2L2', driver)
)

tree_gg_concat <- lapply(names(tree_gg_list_final_2), tree_gg_list_final_2= tree_gg_list_final_2, function(x, tree_gg_list_final_2){
  updt_tree <- tree_gg_list_final_2[[x]]$tree %>% 
    dplyr::mutate(node_color = ifelse(invasive == 'invasive', "11", 
                                      ifelse(node_color == 9, "9", node_color)))
  return(updt_tree)
})

names(tree_gg_concat) <- names(tree_gg_list_final_2)
class(tree_gg_concat) <- "multiPhylo"


## drawing part 
ras <- c("C07","C08","C06","C10") # minimally, widely, adenoma // Excluded C10 ,"A03"
non_ras <- c("C03","C04","C05","C01","C02","C09")  # Excluded ,"A01","A02"
FTA <- c('A01','A03','A02')

# color 
invasive_color = c("invasive" = "white","non_invasive" ="black")
node_color = c("1" = '#FFFF00',"2"= '#00FFFF'  , 
               "3" = '#CCFF00',"4" = '#FFA500', "9"= "white", '11'= 'red' )

whole <- c(non_ras, FTA, ras)

pdf('Fig3a_phylogenetic_trees.pdf',width = 2.8, height = 14)

# 1. build the ggtree object
p <- ggtree(tree_gg_concat[whole], ladderize = FALSE)

patient_order <- c('C09','C05','C04','C03','C02','C01', #  'C01','C02','C03','C04','C05','C06','C09',
                  'A03','A02','A01','C10','C08','C07','C06') #         ,'C07','C08','C10','A01','A02','A03')
# 2. reorder .id in the ggtree data
p$data$.id <- factor(p$data$.id, 
                     levels = patient_order
                     )

# 3. add the remaining ggplot layers
p +
  geom_rootedge(rootedge = unique(tree_gg_concat[whole]$mrca_snv_mut)) +
  geom_point2(aes(subset = (!is.na(label)), fill = node_color), 
              color = 'black', size = 3, shape = 21) +
  geom_text2(aes(x = branch, label = driver), fontface="italic",
             color = 'black', size = 7/.pt, nudge_y = 1.5, angle = 90) +
  geom_text_repel(aes(label = cna_branch),max.overlaps = 12, 
                  color = 'black', size = 7/.pt, nudge_y = 0.5, angle = 90) + # x = x_label_0.6, 
  geom_tiplab(aes(label = ifelse(label2 == "C10-META", "", gsub("^[AC][0-9]{2}","",label2)), color = invasive), 
              size = 6/.pt, offset = -170, angle = 90, hjust = 0.5) +
  theme_tree2() +
  theme(strip.background = element_blank(), 
        strip.text.x     = element_blank(),
        strip.text       = element_blank(),
        panel.spacing    = unit(0, "mm"), 
        legend.position  = "none", 
        axis.text.x      = element_text(size = 5, angle = 90)
        ) +
  scale_color_manual(values = invasive_color) +
  scale_fill_manual(values = node_color) +
  scale_x_continuous(breaks = c(0, 500, 1000, 1500, 2000, 2500),
                     limits = c(0, 3000)) + scale_y_reverse() +
  facet_wrap(~.id,  ncol = 2, dir = "v", scales = 'fixed')

dev.off()


# (FIGURE 3.D) branch-to-trunk ratio by RAS status ------------------------------------

metadata_file = file.path(DATA_DIR, "thyroid_sample_metadata.xlsx")
metadata <- read_excel(metadata_file, sheet = "Analysis_v5", skip= 8)

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

paired_metadata <- metadata %>% dplyr::filter(Tumor == "T" & merged_bam == "F" & High_depth == "T" & Metastasis == "F" ) %>% 
  dplyr::mutate(cancer = ifelse(grepl('FTA',Diagnosis), 'FTA', 
                                ifelse(grepl('FTC',Diagnosis), 'FTC', Diagnosis))) %>% 
  dplyr::select(ID, Diagnosis, cancer, Patient, specific_type, specific_ras, invasive, Tip_snv_mut, Diagnosis_age, Ras, cna_tip_total, Tip_total_mut, MRCA_snv_mut,
                Driver_mut_num, Driver_mut_num_cosmic, Tip_driver_mut_cosmic) %>%
  dplyr::mutate(branch_snv_mut = Tip_snv_mut - MRCA_snv_mut,
                branch_ratio = branch_snv_mut / MRCA_snv_mut)
paired_metadata$Patient <- factor(paired_metadata$Patient, levels = unique(paired_metadata$Patient))

RAS_paired_metadata <- paired_metadata %>%
  group_by(Ras, Patient) %>%   # group by two keys
  mutate(
    branch_ratio_mean = mean(branch_ratio, na.rm = TRUE),
    cna_tip_total_variability = sd(cna_tip_total) / (mean(cna_tip_total) + 1e-6),
    specific_ras = factor(specific_ras,
                          levels = c("F-FTA","F-minimallyFTC","T-minimallyFTC","F-widelyFTC","T-widelyFTC"))
  )  %>% dplyr::filter(!is.na(specific_ras)) %>% 
  ungroup() %>% dplyr::distinct( Ras, Patient, .keep_all = T)

# (an unused summary carried over from 11_fig2.R was removed here; its result was not assigned)




BOXPLOT_WIDTH = 0.6
JITTER_WIDTH = 0.1
JITTER_SIZE = 0.9
FIG3_AXIS_TITLE_SIZE = 7
FIG3_AXIS_TEXT_SIZE=7
FIG3C_annotate_size = 3
FIG3C_Y_MAX = 2
FIG3C_SIGNIF = 2.9
FIG3C_SIGNIF_size = 8
fig3_c <- ggplot(RAS_paired_metadata, aes(x = Ras, y = branch_ratio_mean)) +  
  geom_boxplot(width = BOXPLOT_WIDTH, outlier.shape = NA) +
  geom_jitter(aes(fill = Diagnosis), pch=21, colour = 'black', width = JITTER_WIDTH, size = JITTER_SIZE, alpha = 1, height = 0.08) +
  theme_classic(base_family = 'Arial') +
  scale_x_discrete(
    breaks = c("F","T"),   
    labels = c("RAS\nwt","RAS\nmut")  
  ) +

  # font / axis settings
  theme(
    text = element_text(family = "Arial"),        # Arial font
    axis.text = element_text(size = FIG3_AXIS_TEXT_SIZE),
    axis.title = element_text(size = FIG3_AXIS_TEXT_SIZE)
  ) +
  scale_y_continuous(limits = c(0, FIG3C_Y_MAX), breaks = c(0,1,2), labels = c(0,1,2)) +
  scale_fill_manual(values = c('FTA' = patient_cols[1], 'minimallyFTC' = patient_cols[2], 'widelyFTC' = patient_cols[3])) + 
  ylab('Branch Ratio') + 
  xlab('') +
  
  # comparison brackets and p-values
  geom_signif(comparisons = list(
    c("F", "T")
  ),
  test = "wilcox.test",
 color = 'black', family = 'Arial',
  textsize = FIG3C_SIGNIF_size/.pt,
  )  + # larger p-value font
  theme(legend.position = "none")


ggsave(filename = 'Fig3d_branch_to_trunk_ratio.pdf', path = FIGUREPATH, plot = fig3_c,
       height = 1.7, width = 1.6)



# (FIGURE 3.E) mutation distance vs. physical distance -----------------------------

# mutation distance between sister tips (tips sharing a parent node) in each non-RAS tumor
ras <- c("C07","C08","C06","C10") # RAS-mutant tumors
non_ras <- c("C03","C04","C05","C01","C02","C09")
tree_gg_list_final <- readRDS(file.path(DATA_DIR, "50_additional/sample_mutation_with_drivers/tree_gg_list.rds"))

library(dplyr)
library(data.table)
library(ape)

mut_dist_adjacent <- lapply(non_ras, function(patient) {
  # 1) read the tree and convert to phylo
  ggt <- tree_gg_list_final[[patient]]$tree
  phy <- as.phylo(ggt)
  
  # 2) tip information
  ntip <- length(phy$tip.label)                # tip index: 1..ntip
  tip_ids <- seq_len(ntip)
  tip_labels <- phy$tip.label                  # sample names
  
  # 3) distance matrix between all nodes
  dist_mat <- ape::dist.nodes(phy)
  
  # 4) parent internal node of each tip (edge: [parent, child])
  parent_map <- setNames(phy$edge[,1], phy$edge[,2])   # child -> parent
  tip_parents <- parent_map[as.character(tip_ids)]     # length = number of tips
  
  # 5) group tips that share a parent and build tip pairs within each group
  pairs <- lapply(split(tip_ids, tip_parents), function(tips) {
    if (length(tips) < 2) return(NULL)
    # all tip-tip combinations within the group
    as.data.frame(t(combn(tips, 2)), stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, pairs)
  if (is.null(pairs)) return(NULL)  # NULL if there are no sister tips
  
  colnames(pairs) <- c("i", "j")
  
  # 6) distance for each tip-tip pair
  out <- pairs %>%
    mutate(
      sample1  = tip_labels[i],
      sample2  = tip_labels[j],
      distance = dist_mat[cbind(i, j)]
    ) %>%
    select(sample1, sample2, distance) %>%
    mutate(Patient = patient, .before = 1) %>%
    rename(mut_dist = distance)
  
  out
}) %>% data.table::rbindlist(fill = TRUE)

mut_dist_dat <- mut_dist_adjacent




## CALCULATE THE PHYSICAL DISTANCE FROM SAMPLES 
phy_dist <- readxl::read_excel(metadata_file, sheet = 'distance_mutation')
phy_dist_calc <- function(sample1, sample2, phy_dist){
  # subset the phy_dist
  sam_phy_dist <- phy_dist %>% dplyr::filter(sample %in% c(sample1, sample2)) 
  if(nrow(sam_phy_dist)== 2){
    # extract x difference
    x_diff = diff(sam_phy_dist[['x']])
    # extract y difference 
    y_diff = diff(sam_phy_dist[['y']])
    # extract scale
    scale_1mm <- unique(sam_phy_dist[['scale_5mm']]) / 5
    
    physical_distance <- sqrt(x_diff^2 + y_diff^2) / scale_1mm
  } else{physical_distance <- NA}
  return(physical_distance)
}

mut_phy_dist_dat <- lapply(1:nrow(mut_dist_dat), mut_dist_dat= mut_dist_dat, function(i,mut_dist_dat){
  slice_mut_dist_dat <- mut_dist_dat %>% dplyr::slice(i)
  sample1 = slice_mut_dist_dat[['sample1']]
  sample2 = slice_mut_dist_dat[['sample2']]
  
  slice_mut_dist_dat[['phy_dist_scaled']] <- phy_dist_calc(sample1, sample2, phy_dist)
  return(slice_mut_dist_dat)
}) %>% data.table::rbindlist()


mut_phy_dist_dat <- mut_phy_dist_dat %>% 
  dplyr::mutate(patient = gsub('_.*','', sample1),,
                patient = factor(patient)) %>%
  dplyr::filter(!is.na(phy_dist_scaled)) |> 
  mutate(phy_dist_cm = phy_dist_scaled / 10)

# print the slope (mutation distance per cm) for each patient
for (i_patient in unique(mut_phy_dist_dat$patient)){
  slice_mut <- mut_phy_dist_dat |> filter(patient == i_patient)
  slpe <- lm(mut_dist ~ phy_dist_cm + 0, data = slice_mut)
  print(paste(i_patient, 'slope', slpe$coefficients[["phy_dist_cm"]]))
}


# make class for drawing extended line
lm_left <- function(formula,data,...){
  mod <- lm(formula,data)
  class(mod) <- c('lm_left',class(mod))
  mod
}
lm_left_0 <- function(formula, data, ...) {
  mod <- lm(formula, data)
  class(mod) <- c('lm_left_0', class(mod))  # class name must match exactly
  mod
}

predictdf.lm_left_0 <- function(model, xseq, se, level) {
  init_range <- range(model$model$x)
  xseq <- unique(sort(c(0, xseq[xseq <= init_range[2]])))  # include 0
  ggplot2:::predictdf.default(model, xseq, se, level)
}

distance_plot <- ggplot(mut_phy_dist_dat, aes(x = phy_dist_cm, y = mut_dist)) +
  # confidence band
  stat_smooth(
    data = mut_phy_dist_dat,
    method = "lm_left",
    formula = y ~ 0 + x,
    geom = "ribbon",
    fullrange = TRUE,
    fill = "gray80",
    alpha = 0.5,
    se = TRUE
  ) +
  # line
  stat_smooth(
    data = mut_phy_dist_dat,
    method = "lm_left_0",
    formula = y ~ 0 + x,
    geom = "line",
    fullrange = TRUE,
    col = 'black',
    size = 0.5,
    se = TRUE
  ) +
  geom_smooth(method = lm_left_0, formula = y ~ 0 + x, aes(color = patient), size= 0.5, se = FALSE) +
  geom_point(size = 1.2, aes(shape = patient, color = patient)) + 
  scale_x_continuous(breaks = seq(0,2.5,0.5), labels = seq(0,2.5,0.5), limits = c(0,2.1)) + 
  scale_shape_manual(
    values = c(
      "C03" = 16, # red
      "C01" = 17, # green
      "C05" = 15, # purple
      "C02" = 8, # orange
      "C09" = 7,  # brown
      "C04" = 3
    )
  ) +
  scale_color_manual(
    values = c(
      "C03" = "#E41A1C", # red
      "C04" = "#377EB8", # blue
      "C01" = "#4DAF4A", # green
      "C05" = "#984EA3", # purple
      "C02" = "#FF7F00", # orange
      "C09" = "#A65628"  # brown
    )
  ) +
  theme_classic() +
  ylab("Mutation distance") +
  xlab("Microdissected region distance (cm)") +
  theme(text=element_text(size=7), #change font size of all text
        axis.text=element_text(size=7), #change font size of axis text
        axis.title=element_text(size=7), #change font size of axis titles
        plot.title=element_text(size=7), #change font size of plot title
        legend.text=element_text(size=7), #change font size of legend text
        legend.title=element_text(size=7)) +
  theme(
    legend.position = "bottom",      # bottom centre
    legend.direction = "horizontal"  # legend items in a row
  ) + guides(
    color = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1)
  )

## linear mixed effect model
# fixed effect and random slopes only, no intercept
# fixed-effect slope: 726.6 SNVs per cm
library(lme4)
m_no_int <- lmer(
  mut_dist ~ 0 + phy_dist_cm +
    (0 + phy_dist_cm | patient),
  data = mut_phy_dist_dat
)

summary(m_no_int)

beta <- fixef(m_no_int)["phy_dist_cm"] # 726.64

# dashed line: fixed-effect slope from the mixed model
distance_plot2 <- 
  distance_plot +
  geom_abline(
    intercept = 0,
    slope     = beta,
    linetype  = "dashed",
    color     = "black",
    size      = 0.7
  ) + 
  theme(
    legend.position = "bottom",      # bottom centre
    legend.direction = "horizontal"  # legend items in a row
  )



pdf('Fig3e_mutation_vs_physical_distance.pdf', width = 3.5, height = 2.5)
distance_plot2
dev.off()
