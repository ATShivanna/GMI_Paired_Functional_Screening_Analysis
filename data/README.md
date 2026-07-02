# Data input instructions

Input data files are not included in this repository release because the raw screening datasets are available from the corresponding author upon reasonable request, and representative processed datasets are provided with the manuscript as Supplementary Data S1 and Supplementary Data S2.

To run the full analysis workflow locally, create the following directory structure after cloning the repository:

```text
data_raw/
└── cases/
```

Expected local input files include:

- `data_raw/cases/` — raw plate-reader Excel files for paired PBS and GMI screening conditions.
- `data_raw/drug_lookup.xlsx` — compound and assay-well lookup workbook with sheets used by the scripts, including `Drug_Library`, `Drug_Annotation`, and `Assay_Well_Lookup`.
- `data/Supplementary_Data_S2.xlsx` — workbook used by `scripts/09_regenerate_supp_figures_s9_s10.py` to regenerate Supplementary Figures S9 and S10.

The `data_raw/` directory is intentionally excluded from version control by `.gitignore`.
