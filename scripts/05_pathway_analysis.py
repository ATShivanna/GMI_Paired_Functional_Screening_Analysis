import pandas as pd

# ==============================
# INPUTS (your correct files)
# ==============================
ANNOTATION_FILE = "./data_raw/drug_lookup.xlsx"
PRIMARY_DELTA = "./drug_level_output_all178/drug_matrix_DELTA_all178.csv"
SECONDARY_DELTA = "./drug_level_output_secondary_all178/drug_matrix_DELTA_secondary_all178.csv"

# ==============================
# LOAD
# ==============================
ann = pd.read_excel(ANNOTATION_FILE, sheet_name="Drug_Annotation")

p_delta = pd.read_csv(PRIMARY_DELTA, index_col=0)
s_delta = pd.read_csv(SECONDARY_DELTA, index_col=0)

# ==============================
# MAP drug → pathway
# ==============================
mapping = ann[["Drug_ID", "Primary_Target_or_Pathway"]].dropna()

mapping = mapping.set_index("Drug_ID")

# ==============================
# PRIMARY aggregation
# ==============================
p = p_delta.copy()
p["Pathway"] = p.index.map(mapping["Primary_Target_or_Pathway"])

p = p.dropna(subset=["Pathway"])

p_mean = p.groupby("Pathway").mean(numeric_only=True)
p_mean["Delta"] = p_mean.mean(axis=1)

p_out = p_mean[["Delta"]].reset_index()

# ==============================
# SECONDARY aggregation
# ==============================
s = s_delta.copy()
s["Pathway"] = s.index.map(mapping["Primary_Target_or_Pathway"])

s = s.dropna(subset=["Pathway"])

s_mean = s.groupby("Pathway").mean(numeric_only=True)
s_mean["Delta"] = s_mean.mean(axis=1)

s_out = s_mean[["Delta"]].reset_index()

# ==============================
# SAVE
# ==============================
p_out.to_csv("pathway_output_primary.csv", index=False)
s_out.to_csv("pathway_output_secondary.csv", index=False)

print("✅ Pathway files generated")
