import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# =========================
# GLOBAL STYLE
# =========================
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams["font.size"] = 11
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.labelsize"] = 11
plt.rcParams["ytick.labelsize"] = 6
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# =========================
# PATHS
# =========================
LOOKUP_FILE = "./data_raw/drug_lookup.xlsx"

PRIMARY_DIR = "./drug_level_output_all178"
SECONDARY_DIR = "./drug_level_output_secondary_all178"

OUTDIR = "./heatmap_output_178_final"
os.makedirs(OUTDIR, exist_ok=True)

# =========================
# LOAD DRUG ANNOTATION
# =========================
drug_lib = pd.read_excel(LOOKUP_FILE, sheet_name="Drug_Library")
drug_lib = drug_lib.sort_values("Library_Order")

id_to_short = dict(zip(drug_lib["Drug_ID"].astype(str),
                       drug_lib["Drug_Name_Short"].astype(str)))

id_to_full = dict(zip(drug_lib["Drug_ID"].astype(str),
                      drug_lib["Drug_Name_Full"].astype(str)))

# =========================
# LOAD MATRICES
# =========================
def load_data(base_dir, prefix):
    pbs = pd.read_csv(f"{base_dir}/drug_matrix_PBS_{prefix}.csv", index_col=0)
    gmi = pd.read_csv(f"{base_dir}/drug_matrix_GMI_{prefix}.csv", index_col=0)
    delta = pd.read_csv(f"{base_dir}/drug_matrix_DELTA_{prefix}.csv", index_col=0)
    return pbs, gmi, delta

pbs_p, gmi_p, delta_p = load_data(PRIMARY_DIR, "all178")
pbs_s, gmi_s, delta_s = load_data(SECONDARY_DIR, "secondary_all178")

# =========================
# SORT + LABEL
# =========================
def prepare_matrix(pbs, gmi, delta):
    order = delta.mean(axis=1, skipna=True).sort_values(na_position="last").index

    pbs = pbs.loc[order]
    gmi = gmi.loc[order]
    delta = delta.loc[order]

    labels = [id_to_short.get(x, x) for x in order]

    return pbs, gmi, delta, labels

pbs_p, gmi_p, delta_p, labels_p = prepare_matrix(pbs_p, gmi_p, delta_p)
pbs_s, gmi_s, delta_s, labels_s = prepare_matrix(pbs_s, gmi_s, delta_s)

# =========================
# COLORMAPS (GREY FOR NaN)
# =========================
cmap_viab = plt.cm.viridis.copy()
cmap_viab.set_bad(color='lightgrey')

cmap_delta = plt.cm.bwr_r.copy()
cmap_delta.set_bad(color='lightgrey')

# =========================
# PLOT FUNCTION
# =========================
def plot_heatmap(pbs, gmi, delta, labels, title, filename):

    fig, axes = plt.subplots(
        1, 3,
        figsize=(14, 20),
        gridspec_kw={"width_ratios": [1, 1, 1.05], "wspace": 0.3}
    )

    # Mask NaN → becomes grey
    pbs_mask = np.ma.masked_invalid(pbs.values)
    gmi_mask = np.ma.masked_invalid(gmi.values)
    delta_mask = np.ma.masked_invalid(delta.values)

    vmax = 100
    max_delta = np.nanmax(np.abs(delta.values))
    if np.isnan(max_delta) or max_delta == 0:
        max_delta = 1

    # ----- PBS -----
    im1 = axes[0].imshow(
        pbs_mask,
        aspect='auto',
        cmap=cmap_viab,
        vmin=0,
        vmax=vmax
    )
    axes[0].set_title("A. PBS")
    axes[0].set_ylabel("Drug")
    axes[0].set_xlabel("Case")
    axes[0].set_xticks(range(len(pbs.columns)))
    axes[0].set_xticklabels(pbs.columns, rotation=45)
    axes[0].set_yticks(range(len(labels)))
    axes[0].set_yticklabels(labels)

    # ----- GMI -----
    im2 = axes[1].imshow(
        gmi_mask,
        aspect='auto',
        cmap=cmap_viab,
        vmin=0,
        vmax=vmax
    )
    axes[1].set_title("B. GMI")
    axes[1].set_xlabel("Case")
    axes[1].set_xticks(range(len(gmi.columns)))
    axes[1].set_xticklabels(gmi.columns, rotation=45)
    axes[1].set_yticks([])

    # ----- DELTA -----
    im3 = axes[2].imshow(
        delta_mask,
        aspect='auto',
        cmap=cmap_delta,
        vmin=-max_delta,
        vmax=max_delta
    )
    axes[2].set_title("C. Δ(GMI − PBS)")
    axes[2].set_xlabel("Case")
    axes[2].set_xticks(range(len(delta.columns)))
    axes[2].set_xticklabels(delta.columns, rotation=45)
    axes[2].set_yticks([])

    # Colorbars
    cbar1 = fig.colorbar(im1, ax=axes[0], fraction=0.046, pad=0.02)
    cbar1.set_label("Min viability (%)")

    cbar2 = fig.colorbar(im2, ax=axes[1], fraction=0.046, pad=0.02)
    cbar2.set_label("Min viability (%)")

    cbar3 = fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.02)
    cbar3.set_label("Δ viability (%)")

    # Title
    fig.suptitle(title, fontsize=18, fontweight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.98])

    plt.savefig(f"{OUTDIR}/{filename}.png", dpi=600)
    plt.savefig(f"{OUTDIR}/{filename}.pdf")

    plt.show()

# =========================
# RUN PLOTS
# =========================
plot_heatmap(
    pbs_p, gmi_p, delta_p, labels_p,
    "Primary analysis set (178 compounds)",
    "Heatmap_Primary_178_FINAL"
)

plot_heatmap(
    pbs_s, gmi_s, delta_s, labels_s,
    "Secondary analysis set (178 compounds)",
    "Heatmap_Secondary_178_FINAL"
)

print("\nDONE — heatmaps with grey 'no data' generated.")
