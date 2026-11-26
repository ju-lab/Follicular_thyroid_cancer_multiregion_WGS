# author: taewoo
# date: Mar 28 2025
# description: shared libraries, fonts and colour palettes for the figure scripts

suppressMessages(library(dplyr));
suppressMessages(library(readxl))
suppressMessages(library(patchwork))
suppressMessages(library(ape))
suppressMessages(library(ggtree))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(data.table))
suppressMessages(library(treeio))
suppressMessages(library(openxlsx))
suppressMessages(library(grid))
suppressMessages(library(ggpubr))
library(extrafont)
library(stringr)
library(scales)

# PDF text family ---------------------------------------------------------

# Register Arial for PDF output (font_import writes to the extrafont database, so do it only once)
if (!"Arial" %in% fonts()) font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")

# Sample metadata ---------------------------------------------------------

## Analysis_v5: `ID` is the sample name used for the sequencing data (e.g. C01-1) and `label` is the
## region name used in the other sheets and the data files (e.g. C01_P_9). The figure scripts look
## data up by region name, so `ID` is set to the region name (C01-P-9; for FTA, A01-7_8) and the
## sample name is kept as `sample_id`.
read_sample_metadata <- function(file) {
  readxl::read_excel(file, sheet = "Analysis_v5", skip = 8) %>%
    dplyr::filter(!is.na(ID)) %>%
    dplyr::rename(sample_id = ID) %>%
    dplyr::mutate(ID = ifelse(grepl("^A", label), sub("_", "-", label), gsub("_", "-", label)))
}

# Color palettes ----------------------------------------------------------

## diagnosis colours (FTA, minimally invasive FTC, widely invasive FTC)
patient_cols <- c('#507EBE','#F5C34C', '#A3332D')

## Fig. 2a heatmap: region annotation
region = c('noninvasive' = 'white', 'metastasis' = '#E69F00', 'invasive' = '#E41A1C')

## Fig. 2a heatmap: driver mutation types (colours as in TCGA lung adenocarcinoma, Nature 2014)
missense_hexcode <- rgb(50, 108, 173, maxColorValue = 255) #"#326CAD"
promoter_hexcode<- rgb(231, 136, 35, maxColorValue = 255) #"#E78823"
nonsense_hexcode <- rgb(212, 0, 34, maxColorValue = 255) #'#D40022'
frameshift_hexcode  <- rgb(60, 37, 128, maxColorValue = 255) #'#3C2580'
missense_nonsense_hexcode <- rgb(139, 9, 25, maxColorValue = 255) #'#8B0919'

af_cols = c(frameshift_hexcode, missense_hexcode, nonsense_hexcode, promoter_hexcode, missense_nonsense_hexcode) # frameshift, missense, nonsense, promoter,  missense & nonsense

## Fig. 2a heatmap: copy-number colours
cn_cols = c('blue','white','red')
