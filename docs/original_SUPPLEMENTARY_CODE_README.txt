Supplementary Code – EVETT–GMI Paired Functional Screening Study

This repository contains scripts used for data processing, analysis, and figure generation.

Workflow:

1. 01_extract_normalize.py
   Extracts raw CTG plate data and performs plate-wise normalization.

2. 02_build_drug_matrices.py
   Constructs drug-level matrices from plate-level data.

3. 03_compute_delta.py
   Computes Δ(GMI − PBS) for each drug and case.

4. 04_top_hits_analysis.py
   Identifies top sensitized and attenuated compounds.

5. 05_pathway_analysis.py
   Maps drugs to pathways and aggregates pathway-level responses.

6. 06_clustering.py
   Performs hierarchical clustering across cases based on pathway profiles.

7. 07_statistics.py
   Computes pathway-level statistics (Kruskal–Wallis test).

8. 08_plot_figures.py
   Generates figures used in the manuscript.

Dependencies:
Python 3.9+
pandas
numpy
matplotlib
seaborn
scipy
scikit-learn

All analyses were performed on the final paired dataset (n = 10 cases).
