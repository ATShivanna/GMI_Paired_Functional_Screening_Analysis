import pandas as pd
import numpy as np
import os
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score

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
    df = pd.DataFrame({"Drug_ID": d.index.astype(str).str.strip(), "Delta": d.values})
    df["Pathway"] = df["Drug_ID"].map(drug_to_pathway)
    df = df.dropna(subset=["Delta", "Pathway"]).copy()

    pw = df.groupby("Pathway")["Delta"].mean()
    pw.index = [shorten_pathway(x) for x in pw.index]
    case_pathway[case] = pw

matrix = pd.DataFrame(case_pathway).T.fillna(0)

# Z-score
X = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
X = X.replace([np.inf, -np.inf], 0).fillna(0)

# cluster
dist = pdist(X.values, metric="euclidean")
Z = linkage(dist, method="ward")
cluster_ids = fcluster(Z, t=3, criterion="maxclust")

# silhouette
score = silhouette_score(X.values, cluster_ids, metric="euclidean")

out = pd.DataFrame({
    "Case": X.index,
    "Cluster_ID": cluster_ids
}).sort_values("Case")

csv_out = os.path.join(OUTDIR, "case_cluster_assignments_for_silhouette.csv")
txt_out = os.path.join(OUTDIR, "silhouette_score.txt")

out.to_csv(csv_out, index=False)

with open(txt_out, "w") as f:
    f.write(f"Silhouette score (k=3): {score:.4f}\n")

print(f"Silhouette score (k=3): {score:.4f}")
print("Saved:")
print(csv_out)
print(txt_out)
