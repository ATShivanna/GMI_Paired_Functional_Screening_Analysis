import pandas as pd
import numpy as np
import os
from scipy.stats import kruskal
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

# =============================
# PATHS
# =============================
LOOKUP_FILE = "./data_raw/drug_lookup.xlsx"
PRIMARY_DELTA = "./drug_level_output_all178/drug_matrix_DELTA_all178.csv"
SECONDARY_DELTA = "./drug_level_output_secondary_all178/drug_matrix_DELTA_secondary_all178.csv"
OUTDIR = "./tables"
os.makedirs(OUTDIR, exist_ok=True)

# =============================
# HELPERS
# =============================
def find_col(cols, include_keywords, exclude_keywords=None):
    exclude_keywords = exclude_keywords or []
    cols_l = {c: str(c).strip().lower() for c in cols}
    for c, lc in cols_l.items():
        if all(k in lc for k in include_keywords) and not any(k in lc for k in exclude_keywords):
            return c
    return None

def shorten_pathway(p):
    if pd.isna(p):
        return p
    p = str(p).strip()
    replacements = {
        "SMO / developmental signaling": "SMO/dev",
        "DNA synthesis / replication stress": "DNA replication",
        "Androgen / estrogen / endocrine signaling": "Endocrine",
        "DNA repair / PARP": "DNA repair/PARP",
        "CDK4/6 / cell-cycle control": "CDK4/6-cell cycle",
        "Angiogenesis / proliferation signaling": "Angiogenesis/prolif.",
        "BRAF / MEK signaling": "BRAF/MEK",
        "PI3K / AKT / mTOR axis": "PI3K/AKT/mTOR",
        "Diverse emerging targets": "Emerging targets",
        "DNA damage / alkylation / topoisomerase / broad cytotoxicity": "Broad cytotoxicity",
        "Mitotic spindle / tubulin": "Mitotic spindle",
        "EGFR / ALK / MET / VEGFR / FGFR / BTK / FLT3 / KIT / Src family": "RTK-family",
        "ROS / lipid peroxidation / metabolic stress": "ROS/metabolic",
        "Immune activation / repurposed / metabolic modulation": "Immune/metabolic",
        "Diverse or unresolved target space": "Unresolved targets",
        "HDAC / chromatin regulation": "HDAC/chromatin",
        "Proteasome inhibition": "Proteasome",
    }
    return replacements.get(p, p)

# =============================
# LOAD ANNOTATION
# =============================
ann = pd.read_excel(LOOKUP_FILE, sheet_name="Drug_Annotation")
ann_cols = list(ann.columns)

drug_id_col = (
    find_col(ann_cols, ["drug", "id"]) or
    find_col(ann_cols, ["compound", "id"])
)
pathway_col = (
    find_col(ann_cols, ["primary", "target", "pathway"]) or
    find_col(ann_cols, ["pathway"]) or
    find_col(ann_cols, ["target"])
)

if drug_id_col is None or pathway_col is None:
    raise ValueError(f"Could not identify Drug_ID/pathway columns in Drug_Annotation. Found: {ann_cols}")

ann = ann[[drug_id_col, pathway_col]].copy()
ann.columns = ["Drug_ID", "Pathway"]
ann["Drug_ID"] = ann["Drug_ID"].astype(str).str.strip()
ann["Pathway"] = ann["Pathway"].astype(str).str.strip()
ann = ann.dropna(subset=["Pathway"])
drug_to_pathway = dict(zip(ann["Drug_ID"], ann["Pathway"]))

# =============================
# LOAD DELTA MATRICES
# =============================
p_delta = pd.read_csv(PRIMARY_DELTA, index_col=0)
s_delta = pd.read_csv(SECONDARY_DELTA, index_col=0)
delta_all = pd.concat([p_delta, s_delta], axis=1)

# =============================
# BUILD CASE × PATHWAY MATRIX
# =============================
case_pathway = {}

for case in delta_all.columns:
    d = delta_all[case].copy()
    df = pd.DataFrame({
        "Drug_ID": d.index.astype(str).str.strip(),
        "Delta": d.values
    })
    df["Pathway"] = df["Drug_ID"].map(drug_to_pathway)
    df = df.dropna(subset=["Delta", "Pathway"]).copy()

    pw = df.groupby("Pathway")["Delta"].mean()
    pw.index = [shorten_pathway(x) for x in pw.index]
    case_pathway[case] = pw

matrix = pd.DataFrame(case_pathway).T.fillna(0)

# save case-pathway matrix too
matrix_csv = os.path.join(OUTDIR, "case_by_pathway_matrix_for_stats.csv")
matrix.to_csv(matrix_csv)

# =============================
# Z-SCORE + CLUSTERING
# =============================
X = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
X = X.replace([np.inf, -np.inf], 0).fillna(0)

dist = pdist(X.values, metric="euclidean")
Z = linkage(dist, method="ward")

cluster_ids = fcluster(Z, t=3, criterion="maxclust")
case_to_cluster = {X.index[i]: f"C{cluster_ids[i]}" for i in range(len(X.index))}

# =============================
# LONG FORMAT FOR STATS
# =============================
long_df = matrix.reset_index().rename(columns={"index": "Case"})
long_df = long_df.melt(id_vars="Case", var_name="Pathway", value_name="Delta")
long_df["Cluster"] = long_df["Case"].map(case_to_cluster)

# save cluster assignment
cluster_df = pd.DataFrame({
    "Case": list(case_to_cluster.keys()),
    "Cluster": list(case_to_cluster.values())
}).sort_values("Case")
cluster_csv = os.path.join(OUTDIR, "case_clusters.csv")
cluster_df.to_csv(cluster_csv, index=False)

# =============================
# KRUSKAL-WALLIS STATS
# =============================
rows = []

for pathway, sub in long_df.groupby("Pathway"):
    c1 = sub[sub["Cluster"] == "C1"]["Delta"].dropna()
    c2 = sub[sub["Cluster"] == "C2"]["Delta"].dropna()
    c3 = sub[sub["Cluster"] == "C3"]["Delta"].dropna()

    mean_c1 = c1.mean() if len(c1) > 0 else np.nan
    mean_c2 = c2.mean() if len(c2) > 0 else np.nan
    mean_c3 = c3.mean() if len(c3) > 0 else np.nan

    sd_c1 = c1.std() if len(c1) > 1 else np.nan
    sd_c2 = c2.std() if len(c2) > 1 else np.nan
    sd_c3 = c3.std() if len(c3) > 1 else np.nan

    valid_groups = [g for g in [c1, c2, c3] if len(g) >= 2]
    if len(valid_groups) >= 2:
        stat, pval = kruskal(*valid_groups)
    else:
        pval = np.nan

    if pd.isna(pval):
        sig = "NA"
    elif pval < 0.001:
        sig = "***"
    elif pval < 0.01:
        sig = "**"
    elif pval < 0.05:
        sig = "*"
    else:
        sig = "ns"

    rows.append({
        "Pathway": pathway,
        "Mean_C1": round(mean_c1, 2) if not pd.isna(mean_c1) else np.nan,
        "Mean_C2": round(mean_c2, 2) if not pd.isna(mean_c2) else np.nan,
        "Mean_C3": round(mean_c3, 2) if not pd.isna(mean_c3) else np.nan,
        "SD_C1": round(sd_c1, 2) if not pd.isna(sd_c1) else np.nan,
        "SD_C2": round(sd_c2, 2) if not pd.isna(sd_c2) else np.nan,
        "SD_C3": round(sd_c3, 2) if not pd.isna(sd_c3) else np.nan,
        "Kruskal_p": round(pval, 5) if not pd.isna(pval) else np.nan,
        "Significance": sig
    })

out = pd.DataFrame(rows).sort_values(by="Kruskal_p", ascending=True, na_position="last")

# =============================
# SAVE
# =============================
csv_path = os.path.join(OUTDIR, "pathway_cluster_statistics.csv")
xlsx_path = os.path.join(OUTDIR, "pathway_cluster_statistics.xlsx")

out.to_csv(csv_path, index=False)
out.to_excel(xlsx_path, index=False)

# =============================
# PRINT PREVIEW
# =============================
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
print("\n=== Cluster Pathway Statistics Preview ===\n")
print(out.head(20).to_string(index=False))

print("\nSaved:")
print(csv_path)
print(xlsx_path)
print(cluster_csv)
print(matrix_csv)
