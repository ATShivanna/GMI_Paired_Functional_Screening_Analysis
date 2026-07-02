import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# =============================
# PATHS
# =============================
PRIMARY_DELTA = "./drug_level_output_all178/drug_matrix_DELTA_all178.csv"
SECONDARY_DELTA = "./drug_level_output_secondary_all178/drug_matrix_DELTA_secondary_all178.csv"
OUTDIR = "./figures"
os.makedirs(OUTDIR, exist_ok=True)

# =============================
# STYLE
# =============================
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# =============================
# LOAD
# =============================
d1 = pd.read_csv(PRIMARY_DELTA, index_col=0)
d2 = pd.read_csv(SECONDARY_DELTA, index_col=0)
delta = pd.concat([d1, d2], axis=1)

all_vals = delta.values.flatten()
all_vals = all_vals[~np.isnan(all_vals)]

case_means = delta.mean(axis=0)

# =============================
# SIMPLE KDE (no seaborn needed)
# =============================
def kde_estimate(x, grid, bw=None):
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 2:
        return np.zeros_like(grid)

    if bw is None:
        std = np.std(x, ddof=1)
        iqr = np.subtract(*np.percentile(x, [75, 25]))
        sigma = min(std, iqr / 1.34) if iqr > 0 else std
        bw = 0.9 * sigma * n ** (-1 / 5) if sigma > 0 else 1.0

    if bw <= 0:
        bw = 1.0

    diff = (grid[:, None] - x[None, :]) / bw
    dens = np.exp(-0.5 * diff**2).sum(axis=1) / (n * bw * np.sqrt(2 * np.pi))
    return dens

# =============================
# FIGURE
# =============================
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

# A. pooled distribution
counts, bins, _ = axes[0].hist(
    all_vals,
    bins=50,
    edgecolor="black",
    alpha=0.85
)
axes[0].axvline(0, linestyle="--", linewidth=1)

# KDE scaled to histogram
grid = np.linspace(all_vals.min(), all_vals.max(), 400)
dens = kde_estimate(all_vals, grid)
bin_width = bins[1] - bins[0]
dens_scaled = dens * len(all_vals) * bin_width
axes[0].plot(grid, dens_scaled, linewidth=2)

axes[0].set_title("A. Pooled Δ(GMI−PBS) distribution", fontweight="bold")
axes[0].set_xlabel("Δ(GMI−PBS) viability (%)")
axes[0].set_ylabel("Count")

# B. per-case mean delta
axes[1].bar(range(len(case_means)), case_means.values)
axes[1].axhline(0, linestyle="--", linewidth=1)
axes[1].set_xticks(range(len(case_means)))
axes[1].set_xticklabels(case_means.index, rotation=45, ha="right")
axes[1].set_title("B. Mean Δ(GMI−PBS) by case", fontweight="bold")
axes[1].set_ylabel("Mean Δ(GMI−PBS) viability (%)")

fig.suptitle(
    "Supplementary. Global distribution of GMI-associated response shifts",
    fontsize=15,
    fontweight="bold"
)

plt.tight_layout(rect=[0, 0, 1, 0.96])

pdf_out = os.path.join(OUTDIR, "Supp_delta_distribution_FINAL.pdf")
png_out = os.path.join(OUTDIR, "Supp_delta_distribution_FINAL.png")
plt.savefig(pdf_out, bbox_inches="tight")
plt.savefig(png_out, dpi=600, bbox_inches="tight")
plt.close()

print("Saved:")
print(pdf_out)
print(png_out)
