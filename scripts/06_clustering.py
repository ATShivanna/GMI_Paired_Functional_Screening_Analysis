import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ======================
# INPUT FILES
# ======================
PRIMARY_FILE = "./pathway_output_primary.csv"
SECONDARY_FILE = "./pathway_output_secondary.csv"
OUTDIR = "./figures"
os.makedirs(OUTDIR, exist_ok=True)

# ======================
# LOAD
# ======================
p = pd.read_csv(PRIMARY_FILE)
s = pd.read_csv(SECONDARY_FILE)

p.columns = ["Pathway", "Primary"]
s.columns = ["Pathway", "Secondary"]

# merge
df = pd.merge(p, s, on="Pathway", how="outer")

# ======================
# OPTIONAL: normalize
# ======================
df = df.set_index("Pathway")

# ======================
# CLUSTER HEATMAP
# ======================
plt.figure(figsize=(6, 8))

sns.clustermap(
    df,
    cmap="coolwarm",
    center=0,
    linewidths=0.5,
    col_cluster=True,
    row_cluster=True,
    figsize=(6, 10)
)

plt.savefig("./figures/Supp_case_pathway_clustering.pdf")
plt.savefig("./figures/Supp_case_pathway_clustering.png", dpi=300)

print("✅ Clustering figure saved")
