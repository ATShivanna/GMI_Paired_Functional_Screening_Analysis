import pandas as pd
import numpy as np
import os

# =========================
# PATHS
# =========================
BASE_DIR = "./data_raw"
NORM_DIR = "./normalized_output"
LOOKUP_FILE = "./data_raw/drug_lookup.xlsx"

OUTPUT_DIR = "./drug_level_output_all178"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# LOAD LOOKUP TABLES
# =========================
assay_lookup = pd.read_excel(LOOKUP_FILE, sheet_name="Assay_Well_Lookup")
drug_lib = pd.read_excel(LOOKUP_FILE, sheet_name="Drug_Library")

# Master 178-drug index in library order
drug_lib = drug_lib.sort_values("Library_Order").reset_index(drop=True)
master_drug_ids = drug_lib["Drug_ID"].astype(str).tolist()

# Build mapping: assay well -> Drug_ID
assay_lookup["key"] = assay_lookup["Assay_Row"].astype(str) + assay_lookup["Assay_Col"].astype(str)
well_to_drug = dict(zip(assay_lookup["key"], assay_lookup["Drug_ID"].astype(str)))

# 32-row plate labels
ROW_LABELS = [chr(i) for i in range(ord("A"), ord("Z") + 1)] + ["AA", "AB", "AC", "AD", "AE", "AF"]

# Primary cases
cases = ["A032", "A043", "A156", "A166", "A173", "A174", "A180"]

# =========================
# FUNCTION: MATRIX -> DRUG VECTOR
# =========================
def matrix_to_drug_vector(mat: np.ndarray) -> dict:
    drug_dict = {}

    for i, row_label in enumerate(ROW_LABELS):
        for j in range(48):
            col_num = j + 1
            key = f"{row_label}{col_num}"

            if key not in well_to_drug:
                continue

            drug_id = well_to_drug[key]
            val = mat[i, j]

            if pd.isna(val):
                continue

            drug_dict.setdefault(drug_id, []).append(float(val))

    # Option 2B: MIN viability across wells/doses for each drug
    return {drug_id: np.min(vals) for drug_id, vals in drug_dict.items() if len(vals) > 0}

# =========================
# BUILD FULL 178-DRUG MATRICES
# =========================
pbs_df = pd.DataFrame(index=master_drug_ids, columns=cases, dtype=float)
gmi_df = pd.DataFrame(index=master_drug_ids, columns=cases, dtype=float)

summary_rows = []

for case in cases:
    print(f"Processing {case}")

    pbs_file = f"{NORM_DIR}/{case}_PBS_normalized.csv"
    gmi_file = f"{NORM_DIR}/{case}_GMI_normalized.csv"

    pbs_mat = pd.read_csv(pbs_file, header=None).values
    gmi_mat = pd.read_csv(gmi_file, header=None).values

    pbs_vec = matrix_to_drug_vector(pbs_mat)
    gmi_vec = matrix_to_drug_vector(gmi_mat)

    # fill full matrices
    for drug_id, val in pbs_vec.items():
        pbs_df.loc[drug_id, case] = val
    for drug_id, val in gmi_vec.items():
        gmi_df.loc[drug_id, case] = val

    summary_rows.append({
        "Case": case,
        "PBS_nonmissing": pd.Series(pbs_vec).shape[0],
        "GMI_nonmissing": pd.Series(gmi_vec).shape[0]
    })

# DELTA with NaN preserved
delta_df = gmi_df - pbs_df

# =========================
# ADD OPTIONAL ANNOTATION TABLE
# =========================
annot_cols = ["Drug_ID", "Drug_Name_Short", "Drug_Name_Full", "Drug_Class", "Primary_Target_or_Pathway", "Clinical_Group"]
annot = drug_lib.copy()
for col in annot_cols:
    if col not in annot.columns:
        annot[col] = np.nan
annot = annot[annot_cols]

# =========================
# SAVE
# =========================
pbs_df.to_csv(f"{OUTPUT_DIR}/drug_matrix_PBS_all178.csv")
gmi_df.to_csv(f"{OUTPUT_DIR}/drug_matrix_GMI_all178.csv")
delta_df.to_csv(f"{OUTPUT_DIR}/drug_matrix_DELTA_all178.csv")
annot.to_csv(f"{OUTPUT_DIR}/drug_annotation_all178.csv", index=False)
pd.DataFrame(summary_rows).to_csv(f"{OUTPUT_DIR}/matrix_build_summary_all178.csv", index=False)

# Missingness summary
missing_summary = pd.DataFrame({
    "PBS_missing_n": pbs_df.isna().sum(axis=1),
    "GMI_missing_n": gmi_df.isna().sum(axis=1),
    "DELTA_missing_n": delta_df.isna().sum(axis=1)
}, index=master_drug_ids)
missing_summary.index.name = "Drug_ID"
missing_summary.to_csv(f"{OUTPUT_DIR}/drug_missingness_all178.csv")

print("\nDONE.")
print("Saved to:", OUTPUT_DIR)
print("PBS shape:", pbs_df.shape)
print("GMI shape:", gmi_df.shape)
print("DELTA shape:", delta_df.shape)
print("Total library drugs:", len(master_drug_ids))
print("PBS total missing cells:", int(pbs_df.isna().sum().sum()))
print("GMI total missing cells:", int(gmi_df.isna().sum().sum()))
print("DELTA total missing cells:", int(delta_df.isna().sum().sum()))
