#!/usr/bin/env bash
set -euo pipefail

python scripts/01_extract_normalize.py
python scripts/02_build_drug_matrices.py
python scripts/03_compute_delta.py
python scripts/04_top_hits_analysis.py
python scripts/05_pathway_analysis.py
python scripts/06_clustering.py
python scripts/07_statistics.py
python scripts/07b_silhouette.py
python scripts/08_plot_figures.py
