# Version 1.1.0 — Independent validation release

This release extends the original paired functional screening workflow with reviewer-requested independent validation.

## Added
- Validation analysis for A043, A166, and A180.
- PBS, active GMI, and heat-inactivated GMI conditions.
- Reverse-audited plate QC using Z′ and signal window.
- Complete A–AF assay-well mapping.
- 175-compound original-versus-validation reproducibility analysis.
- 178-compound Active GMI versus heat-inactivated GMI specificity analysis.
- ΔMin and ΔAUC statistics, bootstrap confidence intervals, directional concordance, paired Wilcoxon testing, and BH-FDR correction.
- Supplementary Figures S12–S15 and Supplementary Data S3 source data.

## Compatibility
The original v1.0.x analysis scripts are retained. Validation additions are located in `scripts/validation/`.

## Archival
Zenodo concept DOI: 10.5281/zenodo.21128166

Zenodo will assign a new version-specific DOI when GitHub release v1.1.0 is archived.

## Metadata correction
- Updated the complete nine-author list and affiliations in `CITATION.cff`, `.zenodo.json`, and `README.md`.
- Corresponding authors: Anilkumar T. Shivanna and Ying-Ta Wu.
