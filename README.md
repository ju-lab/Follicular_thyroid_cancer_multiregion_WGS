# Follicular thyroid cancer multiregion WGS

Analysis code for

> Taewoo Kim, Minchul Song, Shinae Yi, Seong Eun Lee, Jin Man Kim, Na Rae Choi, Myeongguk Jeong, Seung-Nam Jung, Zhang Jing, Bon Seok Koo*, Young Seok Ju*, Yea Eun Kang*.
> **Tracing spatial genomic evolution in the clonal history of follicular thyroid neoplasms.**
> *Co-corresponding authors.

## Repository structure

| Script | Figure output (`figures/`) |
|---|---|
| `Thyroid_sources.R` | Shared libraries, fonts and colour palettes (sourced by the figure scripts) |
| `11_fig1.R` | `Fig1b_SNVs_and_signatures.pdf`, `Fig2a_genomic_landscape.pdf` |
| `11_fig2.R` | `Fig2b-e_drivers_and_ploidy.pdf` (exact Wilcoxon rank-sum tests) |
| `11_fig3.R` | `Fig3a_phylogenetic_trees.pdf`, `Fig3d_branch_to_trunk_ratio.pdf`, `Fig3e_mutation_vs_physical_distance.pdf` |
| `11_fig4.R` | `Fig4a_clock_mutation_rate.pdf`, `Fig4c_copy_number_gain_timing.pdf`, `Fig4d_tumor_expansion_time.pdf` |
| `TIMINGANALYSIS_AMPONLY.R` | Timing of copy-number gains relative to point mutations (sourced by `11_fig4.R`) |

Run each script from the repository root, e.g. `Rscript 11_fig1.R`.

## Data (`data/`)

Derived, de-identified data needed to draw the figures. Patient IDs follow the paper (A01–A03, FTA; C01–C10, FTC).

| File | Content |
|---|---|
| `thyroid_sample_metadata.xlsx` | Sample metadata (`Analysis_v5`), point driver mutations (`Driver`), arm-level copy number (`copynumber`), region coordinates on the slide (`distance_mutation`) and copy-number gain timing (`copynumber_timing2`) |
| `50_additional/clean_segments/` | Copy-number segments per sample (Sequenza) |
| `50_additional/sample_mutation_with_drivers/tree_gg_list.rds` | Phylogenetic trees per patient (ggtree objects) |
| `50_additional/timing_output/` | Early SNV rates on amplified segments, per sample |

Raw sequencing data are available as described in the Data availability section of the paper.

## Requirements

R 4.2.3 with Bioconductor 3.16: ComplexHeatmap 2.14.0, ggtree 3.6.2, treeio 1.22.0, ape 5.8.1, circlize 0.4.16, coin 1.4-3, cowplot 1.1.3, data.table 1.17.0, dplyr 1.1.4, extrafont 0.20, GetoptLong 1.0.5, ggplot2 3.5.2, ggpubr 0.6.0, ggrepel 0.9.6, gtable 0.3.6, lme4 1.1-36, openxlsx 4.2.8, patchwork 1.3.0, purrr 1.0.2, readxl 1.4.5, scales 1.3.0, stringr 1.5.1, tidyr 1.3.1. PDF text uses the Arial font. `11_fig4.R` sources `combine_axes()` from [davidearn/plague_growth](https://github.com/davidearn/plague_growth) (GPL-3.0) at a fixed commit, so it needs internet access.
