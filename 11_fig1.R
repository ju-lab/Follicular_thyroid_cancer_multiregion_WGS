# author: taewoo
# date: Oct 12 2025
# description: Fig. 1b (single base substitutions and mutational signatures per sample) and
#              Fig. 2a (genomic landscape heatmap: driver mutations and copy-number alterations)

# Paths: run from the repository root; input data are read from data/ and figures are written to figures/.
CODE_DIR   <- normalizePath(".")
DATA_DIR   <- file.path(CODE_DIR, "data")
FIGUREPATH <- file.path(CODE_DIR, "figures")
dir.create(FIGUREPATH, showWarnings = FALSE)
source(file.path(CODE_DIR, 'Thyroid_sources.R'))
setwd(FIGUREPATH)
# sample metadata organization --------------------------------------------
## March 11 2025
metadata_file = file.path(DATA_DIR, "thyroid_sample_metadata.xlsx")
metadata <- read_excel(metadata_file, sheet = "Analysis_v5", skip= 8)%>% dplyr::filter(!is.na(Patient)) %>%
  dplyr::mutate(Patient = ifelse(grepl('^A',Patient), gsub('A','FTA', Patient), gsub('C','FTC', Patient)))

# patient order (by Diagnosis + Ras)
patient_order <- metadata %>% dplyr::filter(invasive == 'F') %>% 
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



# Samples shown in Fig. 1b and Fig. 2a ------------------------------------------
heatmap_metadata = metadata %>% dplyr::filter( Tumor == "T" & merged_bam == "F" & High_depth == "T") %>% 
  dplyr::mutate(ID = gsub("-","_",ID),
                ID = ifelse(grepl('^A', ID), gsub('A','FTA', ID), 
                                  gsub('C','FTC', ID)),
                Region = ifelse(Metastasis == 'F', ifelse(invasive == 'T', 'invasive', 'noninvasive'), 'metastasis'))

heatmap_metadata <- lapply(unique(metadata$Patient), heatmap_metadata = heatmap_metadata, function(x, heatmap_metadata){
  sub_heatmap_metadata <- heatmap_metadata %>% 
    dplyr::filter( Tumor == "T" & merged_bam == "F" & High_depth == "T" & Patient == x) %>% 
    dplyr::mutate(ID = gsub("-","_",ID), showID = paste0(gsub("_.*","", ID), '-',row_number()))  
  return(sub_heatmap_metadata)
}) %>% data.table::rbindlist()


# Fig. 1b data -------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

# (1) number of single base substitutions per sample (top panel)

df0 <- heatmap_metadata %>%
  mutate(
    showID = dplyr::coalesce(as.character(showID),
                             as.character(label),
                             as.character(Patient),
                             sprintf("S%03d", row_number())),
    Tip_snv_mut = suppressWarnings(as.numeric(Tip_snv_mut)),
    Patient = dplyr::coalesce(as.character(Patient), sub("-.*$", "", showID))
  )

# Sort by Tip_snv_mut (descending) within each patient, keeping the global order
df_plot <- df0 %>%
  group_by(Patient) %>%
  mutate(showID_f = factor(showID, levels = showID)) %>%
  ungroup() %>%
  mutate(
    # strip prefixes, e.g. "FTC08" -> "C08", "FTA01" -> "A01"
    Patient = sub("^FTC", "C", Patient),
    Patient = sub("^FTA", "A", Patient)
  )



# (2) mutational signature proportions per sample (bottom panel)
#     denominator: sum of sig1/sig2/sig5_40/sig13 (set USE_TIP_AS_DENOM = TRUE to divide by Tip_snv_mut)
df <- heatmap_metadata

# fix showID
if (!"showID" %in% names(df)) df$showID <- NA_character_
df <- df %>%
  mutate(
    showID  = dplyr::coalesce(as.character(showID),
                              as.character(label),
                              as.character(Patient),
                              sprintf("S%03d", row_number())),
    # if Patient is missing, use the part of showID before the hyphen as the patient ID
    Patient = dplyr::coalesce(as.character(Patient), sub("-.*$", "", showID))
  )

# (optional) shorten facet labels: FTAxx -> Axx, FTCxx -> Cxx
df <- df %>%
  mutate(Patient = sub("^FTA", "A", Patient),
         Patient = sub("^FTC", "C", Patient))

# coerce to numeric
num_cols <- c("Tip_snv_mut", "sig1", "sig2", "sig5_40", "sig13")
for (cc in num_cols) if (!cc %in% names(df)) df[[cc]] <- NA_real_
df <- df %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))

# --- set showID levels by Tip_snv_mut (descending) within each patient, matching the barplot x order ---
df_bar <- df %>%
  group_by(Patient) %>%
  ungroup()

df <- df %>%
  mutate(showID = factor(showID, levels = df_bar$showID))

# --- signature proportions ---
USE_TIP_AS_DENOM <- FALSE
sig_cols <- c("sig1", "sig2", "sig5_40", "sig13")

df_sig <- df %>%
  mutate(across(all_of(sig_cols), ~ ifelse(is.na(.), 0, .)))

if (USE_TIP_AS_DENOM) {
  df_sig <- df_sig %>%
    mutate(denom = ifelse(is.na(Tip_snv_mut) | Tip_snv_mut == 0, 1, Tip_snv_mut))
} else {
  df_sig <- df_sig %>%
    rowwise() %>%
    mutate(denom = sum(c_across(all_of(sig_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(denom = ifelse(is.na(denom) | denom == 0, 1, denom))
}

df_sig_long <- df_sig %>%
  pivot_longer(cols = all_of(sig_cols), names_to = "Signature", values_to = "Count") %>%
  mutate(Prop = Count / denom) %>%
  mutate(Signature = factor(Signature,
                            levels = c("sig1","sig2","sig5_40","sig13"),
                            labels = c("SBS1","SBS2","SBS5/40","SBS13")))

# ---- 1) per-sample sum of (SBS1 + SBS5/40) and ordering ----
order_tbl <- df_sig_long %>%
  filter(Signature %in% c("SBS1","SBS5/40")) %>%
  group_by(Patient, showID) %>%
  summarise(major_prop = sum(Prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(Patient) %>%
  mutate(showID_ord = factor(showID, levels = showID)) %>%
  ungroup()

# ---- 2) merge SBS2/SBS13 into one group ----
df_major <- df_sig_long %>%
  mutate(Signature2 = case_when(
    Signature %in% c("SBS1","SBS5/40") ~ as.character(Signature),
    TRUE                               ~ "SBS2/13"
  )) %>%
  group_by(Patient, showID, Signature2) %>%
  summarise(Prop = sum(Prop, na.rm = TRUE), .groups = "drop") %>%
  left_join(order_tbl, by = c("Patient","showID")) %>%
  mutate(showID = showID_ord)  # apply the ordering


# ---- 3) palette: two main colours, the rest light grey ----
pal3 <- c("SBS5/40" = "#4A79A8",  # blue
          "SBS1"    = "#5BA94B",  # green
          "SBS2/13"   = "#DDDDDD")  # light grey




# Fig. 1b: bar plot (capped at 2,500) on top of signature proportions ------------

library(grid)  # unit()
library(dplyr)
library(stringr)
library(tidyr)
library(grid)

# shared style
panel_gap <- unit(0.3, "lines")  # slightly smaller gap between panels
y_grid_col <- "grey85"           # colour of the y major grid
axis_lwd    <- 0.6               # axis line width

# top: SNV counts
# sample number shown under each bar, e.g. "A01-1" -> "1"
lab_suffix <- function(x) sub("^[A-Z0-9]+-", "", x)

p_top <- ggplot(df_plot, aes(x = showID_f, y = Tip_snv_mut)) +
  geom_col(width = 0.75) +
  coord_cartesian(ylim = c(0, 2500), clip = "off") +
  labs(x = NULL, y = "Number of single base substitutions") +
  facet_grid(
    cols = vars(Patient),
    scales = "free_x", space = "free_x"
  ) +
  scale_x_discrete(labels = lab_suffix, drop = TRUE) +
  theme_classic(base_size = 8) +
  theme(
    # remove panel border, keep axis lines only
    panel.border       = element_blank(),
    axis.line          = element_line(linewidth = axis_lwd, colour = "black"),
    
    # grid: faint y major lines only, no x grid
    panel.grid.major.y = element_line(colour = y_grid_col, linewidth = 0.25),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    strip.placement    = "outside",
    
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    
    panel.spacing.x    = panel_gap,
    plot.margin        = margin(t = 8, r = 6, b = 2, l = 6)
  )

# bottom: signature proportions
p_bottom <- ggplot(df_major, aes(x = showID, y = Prop, fill = Signature2)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = pal3, breaks = c("SBS1","SBS5/40","SBS2/13")) +
  scale_y_continuous(
    labels  = scales::percent_format(accuracy = 1),
    sec.axis = dup_axis(labels = scales::percent_format(accuracy = 1), name = NULL)
  ) +
  facet_grid(cols = vars(Patient), scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = lab_suffix, drop = TRUE) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  theme_classic(base_size = 8) +
  theme(
    legend.position    = "none",
    strip.background   = element_blank(),
    strip.text.x       = element_blank(),   # hide patient names in the bottom panel
    strip.placement    = "outside",
    
    panel.border       = element_blank(),
    axis.line          = element_line(linewidth = axis_lwd, colour = "black"),
    
    panel.grid.major.y = element_line(colour = y_grid_col, linewidth = 0.25),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    
    axis.text.x        = element_text(angle = 90, vjust = 0.5, size = 7),
    axis.ticks.y.right = element_line(),
    axis.text.y.right  = element_text(),
    
    panel.spacing.x    = panel_gap,
    plot.margin        = margin(t = 2, r = 6, b = 6, l = 6)
  )

# combine (to add a single outer border, use plot.background)
p_combined <- (p_top / p_bottom) + plot_layout(heights = c(2.2, 1), guides = "collect") +
  theme(plot.background = element_rect(fill = NA, colour = NA))  # no outer border
# for a thin outer border, use colour = "black", linewidth = 0.6

ggsave(plot = p_combined,
       filename = "Fig1b_SNVs_and_signatures.pdf",
       path = FIGUREPATH, width = 7.26, height = 5)


# Fig. 2a: genomic landscape heatmap -----------------------------------------------



## METADATA
hat_mat <- heatmap_metadata %>% as.matrix(.)
rownames(hat_mat) = heatmap_metadata$showID 

## color palettes for each metadata 
heatmap_annotation_cols <- list(
  region = region
)

# drawing the color for signature 
sig_bar <- as.data.frame(hat_mat[,'Tip_total_mut']) %>% dplyr::mutate_if(is.character, as.numeric) %>% as.matrix(.) 
# the last sample (C10-META, 35,972 mutations) is far above the axis limit (3,000);
# its bar is capped for display and the actual value is labelled on the figure
sig_bar[nrow(sig_bar), 1]<- 3363


# information of patient
patient_info <- heatmap_metadata %>% dplyr::distinct(Patient, .keep_all = T) %>% 
  tidyr::unite(col =info,  c('Patient','Sex', 'Diagnosis_age'), sep = "\n") %>% dplyr::pull(info) %>% 
  gsub('$','',.) %>% gsub("\nXX\n","\nF,",.) %>% gsub("\nXY\n","\nM,",.) %>% gsub(',','/',.)

sample_info <- heatmap_metadata  %>% 
  tidyr::unite(col =info,  c('Patient','Sex', 'Diagnosis_age'), sep = "\n") %>% dplyr::pull(info)  %>% 
  gsub('$','',.) %>%  gsub("\nXX\n","\nF,",.) %>% gsub("\nXY\n","\nM,",.) %>% gsub(',','/',.)


# DRIVER PART  ------------------------------------------------------------

## Oncoprint 
### driver mutation part 
drivers <- read_excel(metadata_file, sheet = "Driver") 
colnames(drivers) <- colnames(drivers) %>% gsub( "-", "_", .) 
colnames(drivers) <- ifelse(grepl('^A', colnames(drivers)), colnames(drivers) %>% gsub('A','FTA',.), colnames(drivers) %>% gsub('C', 'FTC', .))
driver_mutations <- drivers[["...1"]]
drivers <- drivers[,-1] %>% 
  dplyr::select(heatmap_metadata[["ID"]]) 
colnames(drivers) <- heatmap_metadata$showID
drivers_mat <- as.matrix(drivers)
rownames(drivers_mat) <- driver_mutations

driver_num <- c(c(3,rep(1,6)),c(5,2,2,2,2)) #,rep(1,20)

## sort genes for thyroid genes
oncoprint_genes <- c("NRAS", "KRAS", "HRAS", "TERT", "TSHR", "DICER1", "EZH1",
                     "TFG", "RBM10","BCL3","FAT1","MED12") #,"RNF213","NFE2L2","RPL5","AFF4","MAP3K1","PDE4DIP","A1CF","HIP","FBLN2","PTEN","SETD1B","RECQL4","MUC16","HNF1A","SMARCA4","CDH1","BRD4","CDX2","SRGAP3","TRRAP")

names(af_cols) = c(
  'Frame_Shift_Mutation',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Promoter_Mutation',
  'Nonsense_and_Missense_Mutation'
)

alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  missense_variant = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill =af_cols['Missense_Mutation'], col = NA)),
  promoter_variant = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Promoter_Mutation'], col = NA)),
  frameshift_variant = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Frame_Shift_Mutation'], col = NA)),
  stop_gained = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Nonsense_Mutation'], col = NA)),
  `missense_variant&stop_gained` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Nonsense_and_Missense_Mutation'], col = NA)),
  `missense_variant&splice_region_variant` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Missense_Mutation'], col = NA)),
  `frameshift_variant&splice_region_variant` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Frame_Shift_Mutation'], col = NA)),
  `frameshift_variant&missense_variant` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = af_cols['Frame_Shift_Mutation'], col = NA))
)

## top annotation 
hat = HeatmapAnnotation(
  empty = anno_empty(border = FALSE), 
  info = anno_block(show_name = T, labels = patient_info,  
                    panel_fun = function(index, nm) {
                      grid.rect(gp = gpar(col = "white"))
                      grid.text(nm, 0.5, 0.5, gp = gpar(fontsize = 6))
                    }),
  region = hat_mat[,"Region"],
  mutations = anno_barplot(
    sig_bar,
    border = FALSE,       # no border around anno_barplot
    bar_width = 0.9,      # bar width
    ylim = c(0,3000),
    gp = gpar(fill = "darkgrey", col = NA),  # bar colour: grey, no border
    axis_param = list(
      side = "left",
      gp = gpar(fontsize = 5),
      at = c(0, 1000, 2000, 3000),
      labels = c(0, 1000, 2000, 3000)
    )
  ), 
  col = heatmap_annotation_cols, 
  show_legend = T, 
  annotation_name_side = 'left', 
  annotation_name_gp = gpar(fontsize = 6),
  annotation_name_offset = c(mutations = "1cm"),
  annotation_label = c(
    empty = "",
    info = "Patient\nSex/Age",
    region = "Region",
    mutations = "Number of\nmutations"
  ), 
  border = TRUE, 
  annotation_legend_param = list(
    region = list(title = "Region",
                    # lay out legend items horizontally
                    at = c('noninvasive','invasive', 'metastasis'),
                    border = 'black',
                    labels = c('Noninvasive region','Invasive region', 'Lung metastasis'),
                    labels_gp = gpar(fontsize = 5),
                    grid_width = unit(1.6,"mm"),
                    title_gp = gpar(fontsize = 5)#,grid_width = unit(1.6,"mm")
    )
  ), 
  gap = unit(1.5,"mm") ,
  annotation_height = unit(c(6, 6, 2, 8), rep("mm", 4))
)


# draw the empty annotation together with the existing oncoPrint object
driver_heat <- oncoPrint(
  drivers_mat[oncoprint_genes,],
  name = "driver",
  top_annotation = hat,
  alter_fun = alter_fun,
  show_pct = FALSE,
  alter_fun_is_vectorized = FALSE,
  show_column_names = TRUE,
  row_order = oncoprint_genes,
  row_split = factor(c(rep("RAS", 3), rep("Thyroid driver genes", 4),
                       rep("COSMIC cancer genes", length(oncoprint_genes) - 7)),
                     levels = c('RAS',"Thyroid driver genes", "COSMIC cancer genes")),
  row_title_gp = gpar(fontsize = 6),  # group label font size
  row_title_side = "left",
  row_names_gp = gpar(fontsize = 5, fontface = 'italic'),
  row_gap = unit(0.7, "mm"),
  column_split = factor(sample_info, levels = patient_info),
  column_order = heatmap_metadata$showID,
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = rep(0, length(patient_info))),
  column_gap = unit(1.5, "mm"),
  row_names_side = 'left', 
  right_annotation = rowAnnotation(
    annotation_label = '', 
    count = anno_barplot(
      driver_num, 
      border = FALSE,
      bar_width = 0.9,
      ylim = c(0,6),
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 5),
      width = unit(0.7, "cm"),
      gp = gpar(fill = 1),
      axis_param = list(
        side = "top", 
        labels_rot = 0, 
        gp = gpar(fontsize = 5), 
        at = c(0,4,8), 
        labels = c(0,4,8)
      )
    )
  ), 
  height = unit(4,'cm'), 
  heatmap_legend_param = list(
    title = "Variant", 
    grid_width = unit(2, "mm"),
    border = 'black',
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 5),
    at = c('missense_variant','promoter_variant','frameshift_variant','stop_gained','missense_variant&stop_gained'),
    labels = c('Missense','Promoter','Frameshift','Nonsense','Multiple (missense and nonsense)')
  )
)

# CNV HEATMAP PART --------------------------------------------------------

### CNV part 
copynumber <- read_excel(metadata_file, sheet = "copynumber", skip= 1)
colnames(copynumber) <- gsub("-","_",colnames(copynumber)) %>% gsub('^F','C',.)
colnames(copynumber) <- ifelse(grepl('^A',colnames(copynumber)), gsub('^A','FTA', colnames(copynumber)), 
                               gsub('^C','FTC', colnames(copynumber)))
cnv_num_amp <- c(3,5,2,2,3,3,3,3, 5,4,2,2,4,4,2,1,4,4, 2,2,3,3,5,5,3,6,3,
                 4,4,4,5,2,3,3,3,4,4,4,4,3)
cnv_num_del <- c(0,0,0,0,1,0,0,0, 0,1,0,0,0,0,0,0,0,0,1,0,
                 0,1,0,0,1,0,1,0,0,1,0,2,1,2,0,0,0,0,3,1)
cnv_num <- rbind(cnv_num_amp, cnv_num_del)

#### Matrix build
cnv_mat <- lapply(1:nrow(heatmap_metadata),
                  heatmap_metadata = heatmap_metadata, 
                  copynumber = copynumber,
                  function(i, heatmap_metadata, copynumber){
                    id  <- heatmap_metadata %>% dplyr::slice(i) %>% dplyr::pull(ID) 
                    vec <- lapply(1:nrow(copynumber), function(x, chrs = copynumber[[id]]){
                      str_list <- str_split(chrs[x], pattern = ";") 
                      return( as.numeric(str_list[[1]][1]) / 2)
                    }) %>% unlist()
                    return(vec)
                  }) %>% do.call( cbind, .)
rownames(cnv_mat) <- copynumber[[1]]
colnames(cnv_mat) <- heatmap_metadata$showID

#### Set colors for heatmap
col_fun = colorRamp2(c(0, 1, 2), cn_cols)



#### CNV Heatmap drawing of 
cnv_heat <- Heatmap(cnv_mat, 
                    name = "cnv", 
                    cluster_columns = F, 
                    col = col_fun, 
                    cell_fun = function(j, i, x, y, w, h, fill) {
                      # value from cnv_mat (assumed numeric)
                      value <- cnv_mat[i, j]
                      # 1) draw the full cell border
                      grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "grey"))
                      # 2) fill an inner rectangle at 0.9x size to leave a margin (border)
                      grid.rect(x, y, w * 0.9, h * 0.9, gp = gpar(fill = col_fun(value), col = NA))
                    },
                    show_heatmap_legend = T, #rect_gp = gpar(col = "grey", lwd = 0.9),
                    cluster_rows = F, 
                    row_split = factor(c(rep(1:12, each = 2), 13, 14, 15, rep(16:20, each = 2 ), 21, 22, "X"), levels = c(1:22,"X")), #factor(rownames(cnv_mat), levels = rownames(cnv_mat)), 
                    row_names_side = 'left',
                    row_title = "Copy number alteration",
                    row_title_gp = gpar(fontsize = 6), 
                    row_title_rot = 90,
                    row_names_gp = gpar(fontsize = 5),
                    row_gap = unit(0.5,"mm"),
                    column_split = heatmap_metadata$Patient, 
                    column_title = NULL, 
                    column_order = heatmap_metadata$showID,
                    column_labels = colnames(cnv_mat) %>% gsub('.*-','',.),
                    column_names_rot = 90,
                    column_names_gp = gpar(fontsize = 5),
                    column_gap = unit (1.5, "mm"),
                    heatmap_legend_param =list(
                      # lay out legend items horizontally
                      legend_height = unit(1, "cm"), # legend height
                      title = "Copy number ratio", at = c(0, 0.5, 1, 1.5, 2), 
                      labels = c("0", "0.5", "1", "1.5", "2"),
                      labels_gp = gpar(fontsize = 5), 
                      title_gp = gpar(fontsize = 5),
                      grid_width = unit(2, "mm")
                    ), 
                    height = unit(6.5, 'cm'), 
                    right_annotation = rowAnnotation(
                      annotation_label = '',
                      width = unit(0.7, "cm"),
                      count = anno_barplot(colSums(cnv_num), bar_width = 0.9, border = FALSE, ylim = c(0,6), add_numbers = TRUE,
                                           numbers_gp = gpar(fontsize = 5),  gp = gpar(fill = 1),
                                                         axis_param = list(side = "bottom", labels_rot = 0, 
                                                                           gp = gpar(fontsize = 5), 
                                                                           at = c(0,4,8),
                                                                           labels = c(0,4,8))))
) 


## function for adding block 
library(GetoptLong)  # for the function qq()
group_block_diagnosis = function(group, empty_anno, gp = gpar(), 
                                 label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}

group_block_invasive = function(group, num_noninvasive, num_total, gp = gpar(), 
                                label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_region_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  ## cnv part
  seekViewport(qq("cnv_heatmap_body_1_@{max(group)}"))
  loc3 = deviceLoc(x = unit(0, "npc"), y = unit(1, "npc"))
  seekViewport(qq("cnv_heatmap_body_23_@{max(group)}"))
  loc4 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport("global")
  grid.rect(loc1$x + (loc2$x - loc1$x) * (num_noninvasive / num_total), loc4$y, width = (loc2$x - loc1$x ) * (1-(num_noninvasive / num_total)), height = loc3$y - loc4$y, 
            just = c("left", "bottom"), gp = gp)
  
  ## driver part 
  seekViewport(qq("driver_heatmap_body_1_@{max(group)}"))
  loc5 = deviceLoc(x = unit(0, "npc"), y = unit(1, "npc"))
  seekViewport(qq("driver_heatmap_body_2_@{max(group)}"))
  loc6 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("driver_heatmap_body_3_@{max(group)}"))
  loc7 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport("global")
  grid.rect(loc1$x + (loc2$x - loc1$x) * (num_noninvasive / num_total), loc7$y, width = (loc2$x - loc1$x ) * (1-(num_noninvasive / num_total)), height = loc5$y - loc7$y, 
            just = c("left", "bottom"), gp = gp)
  
  ## mutations part
  seekViewport(qq("annotation_mutations_@{max(group)}"))
  loc7 = deviceLoc(x = unit(0, "npc"), y = unit(1, "npc"))
  seekViewport(qq("annotation_mutations_@{max(group)}"))
  loc8 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport("global")
  grid.rect(loc1$x + (loc2$x - loc1$x) * (num_noninvasive / num_total), loc8$y, width = (loc2$x - loc1$x ) * (1-(num_noninvasive / num_total)), height = loc7$y - loc8$y, 
            just = c("left", "bottom"), gp = gp)
  
  ### metastasis part
  #seekViewport(qq("annotation_Lung_Metastasis_@{max(group)}"))
  #loc9 = deviceLoc(x = unit(0, "npc"), y = unit(1, "npc"))
  #seekViewport(qq("annotation_Lung_Metastasis_@{max(group)}"))
  #loc10 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  #seekViewport("global")
  #grid.rect(loc1$x + (loc2$x - loc1$x) * (num_noninvasive / num_total), loc10$y, width = (loc2$x - loc1$x ) * (1-(num_noninvasive / num_total)), height = loc9$y - loc10$y, 
  #          just = c("left", "bottom"), gp = gp)
}


# DRAW COMBINATED HEATMAP -------------------------------------------------
# 1) set the global font to Arial
# supported options are set with ht_global_opt()
ht_opt(
  heatmap_row_names_gp       = gpar(fontfamily = "Arial"),
  heatmap_column_names_gp    = gpar(fontfamily = "Arial"),
  heatmap_row_title_gp       = gpar(fontfamily = "Arial"),
  heatmap_column_title_gp    = gpar(fontfamily = "Arial"),
  legend_labels_gp   = gpar(fontfamily = "Arial"),
  legend_title_gp    = gpar(fontfamily = "Arial"),
  legend_border              = "black",
  heatmap_border             = FALSE ,
  annotation_border          = 'white'
)

pdf('Fig2a_genomic_landscape.pdf', width = 8)
# concatenate the matrix 
combined_heatmap <-   driver_heat %v% cnv_heat

# draw the combined heatmap with all legends (heatmap and annotation legends) at the bottom
draw(
  combined_heatmap,
  merge_legend = TRUE 
)

# diagnosis annotation
group_block_diagnosis(1:3, "empty", label = "FTA", gp = gpar(fill = adjustcolor(patient_cols[1], alpha.f =0.5), col = NA)) # , 
group_block_diagnosis(4:11, "empty",  label = "Minimally FTC", gp = gpar(fill = adjustcolor(patient_cols[2], alpha.f =0.5), col = NA)) 
group_block_diagnosis(12:13, "empty",  label = "Widely FTC",gp = gpar(fill = adjustcolor(patient_cols[3], alpha.f =0.5), col = NA)) 

# shade invasive regions: per patient (column block), number of non-invasive samples and total samples
group_block_invasive(group = 4,num_noninvasive = 5, num_total = 7,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 5,num_noninvasive = 6, num_total = 10,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 6,num_noninvasive = 8, num_total = 10,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 7,num_noninvasive = 4, num_total = 6,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 8,num_noninvasive = 3, num_total = 4,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 9,num_noninvasive = 2, num_total = 3,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 10,num_noninvasive = 8, num_total = 9,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 11,num_noninvasive = 6, num_total = 7,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 12,num_noninvasive = 7, num_total = 10,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))
group_block_invasive(group = 13,num_noninvasive = 6, num_total = 13,  gp = gpar(fill = adjustcolor("red",alpha.f = 0.1), col = NA))

dev.off()
