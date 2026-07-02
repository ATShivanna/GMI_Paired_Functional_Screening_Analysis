import pandas as pd
import numpy as np
import os

# =========================
# PATHS
# =========================
INPUT_DIR = "./drug_level_output_all178"
LOOKUP_FILE = "./data_raw/drug_lookup.xlsx"

OUTDIR = "./top_hits_output"
os.makedirs(OUTDIR, exist_ok=True)

# =========================
# LOAD DATA
# =========================
pbs = pd.read_csv(f"{INPUT_DIR}/drug_matrix_PBS_all178.csv", index_col=0)
gmi = pd.read_csv(f"{INPUT_DIR}/drug_matrix_GMI_all178.csv", index_col=0)
delta = pd.read_csv(f"{INPUT_DIR}/drug_matrix_DELTA_all178.csv", index_col=0)

# =========================
# LOAD ANNOTATION
# =========================
drug_lib = pd.read_excel(LOOKUP_FILE, sheet_name="Drug_Library")

id_to_short = dict(zip(drug_lib["Drug_ID"].astype(str),
                       drug_lib["Drug_Name_Short"].astype(str)))

id_to_full = dict(zip(drug_lib["Drug_ID"].astype(str),
                      drug_lib["Drug_Name_Full"].astype(str)))

# =========================
# PARAMETERS
# =========================
TOP_N = 10   # change to 20 if needed

# =========================
# PROCESS EACH CASE
# =========================
for case in delta.columns:

    print(f"\nProcessing {case}")

    df = pd.DataFrame({
        "Drug_ID": delta.index,
        "Delta": delta[case],
        "PBS": pbs[case],
        "GMI": gmi[case]
    })

    # Remove NaN (non-evaluable drugs)
    df = df.dropna()

    # Add names
    df["Drug_Short"] = df["Drug_ID"].map(id_to_short)
    df["Drug_Full"] = df["Drug_ID"].map(id_to_full)

    # =========================
    # TOP SENSITIZED (most negative Δ)
    # =========================
    top_sens = df.sort_values("Delta", ascending=True).head(TOP_N)
    top_sens["Direction"] = "Sensitized"

    # =========================
    # TOP ATTENUATED (most positive Δ)
    # =========================
    top_att = df.sort_values("Delta", ascending=False).head(TOP_N)
    top_att["Direction"] = "Attenuated"

    # Combine
    result = pd.concat([top_sens, top_att])

    # Save
    result.to_csv(f"{OUTDIR}/{case}_top_hits.csv", index=False)

    print(f"Saved: {case}_top_hits.csv")

print("\nDONE — Top hits extracted.")

