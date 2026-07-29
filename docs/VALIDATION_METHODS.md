# Independent validation analysis methods

## Design
A043, A166, and A180 were re-screened under PBS, active GMI, and heat-inactivated GMI.

## Plate QC
The audited original workflow used plate-specific high-response controls and the lowest 5% of valid
raw plate values as the low-response reference. Z′ and signal window were calculated from reference
means and standard deviations.

## Normalization
Normalized viability (%) = 100 × (raw − low-reference median) /
(high-reference median − low-reference median), constrained to 0–100%.

## Compound mapping
The complete assay-well lookup includes all 32 rows (A–AF) and 178 compounds.

## Reproducibility population
Three positions contained different compound identities in the validation library:
D024, D060, and D142. They were excluded from direct original-versus-validation comparisons,
leaving 175 identity-matched compounds.

## Endpoints
ΔMin = minimum viability under GMI minus minimum viability under PBS.
ΔAUC = trapezoidal integrated viability under GMI minus PBS across the ordered eight-dose series.

## Statistics
Spearman rank correlations and 3,000-resample bootstrap 95% confidence intervals were used for
reproducibility. ΔAUC directional concordance was also reported. Active GMI versus HI-GMI was
tested within each case using paired Wilcoxon signed-rank tests; the three P values were adjusted
using Benjamini-Hochberg FDR.
