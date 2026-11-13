# Follicular thyroid cancer multiregion WGS

Analysis code for

> Taewoo Kim, Minchul Song, Shinae Yi, Seong Eun Lee, Jin Man Kim, Na Rae Choi, Myeongguk Jeong, Seung-Nam Jung, Zhang Jing, Bon Seok Koo*, Young Seok Ju*, Yea Eun Kang*.
> **Tracing spatial genomic evolution in the clonal history of follicular thyroid neoplasms.**

## Repository structure

| Script | Content |
|---|---|
| `Thyroid_sources.R` | Shared libraries, colour palettes and helper functions (sourced by the figure scripts) |
| `11_fig1.R` | Fig. 1b,c (mutation burden, mutational signatures); genomic landscape heatmap (Fig. 2a) |
| `11_fig2.R` | Fig. 2b–e (point driver mutations and tumor ploidy; exact Wilcoxon rank-sum tests) |
| `11_fig3.R` | Fig. 3 (phylogenetic trees, branch-to-trunk ratio, mutation distance vs. physical distance) |
| `11_fig4.R` | Fig. 4 (clock-like mutation rate, timing of copy-number gains, tumor expansion time) |
| `TIMINGANALYSIS_AMPONLY.R` | Timing of copy-number gains relative to point mutations (sourced by `11_fig4.R`) |
