# Changelog

## v1.1.0

### Independent validation
- Added independent validation analysis for A043, A166, and A180.
- Added PBS, active GMI, and heat-inactivated GMI comparisons.
- Added reverse-audited Z′ and signal-window QC.
- Corrected complete compound mapping across assay rows A–AF.
- Added reproducibility analysis for 175 identity-matched compounds.
- Added within-validation specificity analysis across all 178 compounds.
- Added ΔMin and ΔAUC correlations, bootstrap 95% confidence intervals, directional concordance, paired Wilcoxon tests, and BH-FDR correction.
- Added Supplementary Figures S12–S15 and Supplementary Data S3 source files.
- Standardized Zenodo concept DOI to 10.5281/zenodo.21128166.

## v1.0.4

### Documentation
- Removed broken README banner image reference.
- Updated software badges for the Zenodo-ready release.
- Standardized Python version labeling to Python 3.11.
- Improved citation, data availability, and repository maintainer wording.

### Release
- Repository prepared for Zenodo DOI archiving.
- Repository prepared for the Scientific Reports code availability statement.
- No computational methods or analysis scripts were modified.

All notable changes to this repository are documented here.

## v1.0.1 — Initial public release

Initial GitHub/Zenodo-ready release accompanying the manuscript under revision.

Included:

- Complete Python analysis workflow for paired PBS and GMI functional screening data.
- Scripts for data extraction, viability normalization, Δ(GMI − PBS) calculation, compound prioritization, pathway aggregation, clustering, statistical testing, and figure generation.
- ΔMin–ΔAUC sensitivity analysis figure-regeneration script.
- Regenerated Supplementary Figure S9 and Supplementary Figure S10 outputs.
- Repository metadata files, including `CITATION.cff`, `.zenodo.json`, `LICENSE`, `requirements.txt`, and release notes.
- Clean folder structure without empty placeholder directories.

Notes:

- Raw screening data are not included in this repository release.
- Representative processed datasets are provided with the manuscript as Supplementary Data S1 and Supplementary Data S2.
- Additional datasets are available from the corresponding author upon reasonable request.
