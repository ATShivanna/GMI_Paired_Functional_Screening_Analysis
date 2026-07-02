# Workflow Notes

This document summarizes the computational workflow included in this repository.

## Main analysis pipeline

1. `01_extract_normalize.py` extracts raw plate-reader matrices and performs viability normalization.
2. `02_build_drug_matrices.py` generates compound-level PBS, GMI, and Δ response matrices.
3. `03_compute_delta.py` summarizes Δ(GMI − PBS) response behavior.
4. `04_top_hits_analysis.py` identifies top sensitization- and attenuation-associated compounds.
5. `05_pathway_analysis.py` aggregates compound-level responses by curated pathway annotations.
6. `06_clustering.py` performs hierarchical clustering of pathway-response profiles.
7. `07_statistics.py` performs Kruskal–Wallis testing and BH-FDR correction.
8. `07b_silhouette.py` evaluates supporting cluster structure.
9. `08_plot_figures.py` generates main and supporting figures.

## ΔMin–ΔAUC sensitivity analysis

`09_regenerate_supp_figures_s9_s10.py` regenerates Supplementary Figures S9 and S10 from `Supplementary_Data_S2.xlsx`.

## Data policy

Raw screening files are not included in this public repository release. Representative processed datasets are provided as manuscript supplementary data files. Additional datasets may be requested from the corresponding author.

## Graphical overview

![Workflow overview](workflow_overview.svg)
