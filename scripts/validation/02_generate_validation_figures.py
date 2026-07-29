#!/usr/bin/env python3
"""Generate Supplementary Figures S13-S15 from processed validation data."""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon

EXCLUDED_IDS = {"D024", "D060", "D142"}
CASES = ("A043", "A166", "A180")


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-data", required=True)
    parser.add_argument("--auc-data", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    min_rows = read_csv(args.min_data)
    auc_rows = read_csv(args.auc_data)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    for metric, rows, xcol, ycol, filename in [
        ("ΔMin", min_rows, "Original_Delta", "Validation_Delta_Active_minus_PBS",
         "Supplementary_Figure_S13_DeltaMin_Reproducibility.png"),
        ("ΔAUC", auc_rows, "Original_DeltaAUC", "Validation_DeltaAUC",
         "Supplementary_Figure_S14_DeltaAUC_Reproducibility.png"),
    ]:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
        for ax, case in zip(axes, CASES):
            case_rows = [r for r in rows if r["Case"] == case and r["Drug_ID"] not in EXCLUDED_IDS]
            x = np.asarray([float(r[xcol]) for r in case_rows])
            y = np.asarray([float(r[ycol]) for r in case_rows])
            rho, p = spearmanr(x, y)
            ax.scatter(x, y, s=18, alpha=0.55)
            lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
            ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1)
            ax.axhline(0, linewidth=0.7)
            ax.axvline(0, linewidth=0.7)
            ax.set_title(case)
            ax.set_xlabel(f"Original {metric}")
            if ax is axes[0]:
                ax.set_ylabel(f"Validation {metric}")
            ax.text(0.04, 0.96, f"n = {len(case_rows)}\nρ = {rho:.3f}\nP = {p:.2g}",
                    transform=ax.transAxes, va="top", fontsize=8)
        fig.tight_layout()
        fig.savefig(out / filename, dpi=300, bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    rng = np.random.default_rng(42)
    for ax, case in zip(axes, CASES):
        rows = [r for r in min_rows if r["Case"] == case]
        active = np.asarray([float(r["Validation_Delta_Active_minus_PBS"]) for r in rows])
        hi = np.asarray([float(r["Validation_Delta_HI_minus_PBS"]) for r in rows])
        ax.scatter(rng.normal(1, 0.045, len(active)), active, s=10, alpha=0.25)
        ax.scatter(rng.normal(2, 0.045, len(hi)), hi, s=10, alpha=0.25)
        ax.boxplot([active, hi], positions=[1, 2], widths=0.45, showfliers=False)
        ax.axhline(0, linewidth=0.8)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Active GMI", "HI-GMI"])
        ax.set_title(case)
        if ax is axes[0]:
            ax.set_ylabel("ΔMin versus PBS")
    fig.tight_layout()
    fig.savefig(out / "Supplementary_Figure_S15_Active_vs_HI_GMI.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
