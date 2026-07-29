#!/usr/bin/env python3
"""Recalculate final independent-validation statistics from processed source data.

Reproducibility uses 175 identity-matched compounds.
Within-validation Active GMI versus HI-GMI comparisons use all 178 compounds.
"""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr, wilcoxon

EXCLUDED_IDS = {"D024", "D060", "D142"}
CASES = ("A043", "A166", "A180")


def bootstrap_spearman(x, y, n_boot=3000, seed=42):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(x), len(x))
        xb, yb = x[idx], y[idx]
        if np.unique(xb).size < 2 or np.unique(yb).size < 2:
            continue
        values.append(spearmanr(xb, yb).statistic)
    return np.percentile(values, [2.5, 97.5])


def bh_fdr(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(len(p))
    previous = 1.0
    for i in range(len(p) - 1, -1, -1):
        value = min(previous, ranked[i] * len(p) / (i + 1))
        adjusted[i] = value
        previous = value
    out = np.empty(len(p))
    out[order] = adjusted
    return out


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-data", required=True)
    parser.add_argument("--auc-data", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    min_rows = read_csv(args.min_data)
    auc_rows = read_csv(args.auc_data)
    summary = []
    specificity_p = []
    specificity = []

    for case_index, case in enumerate(CASES):
        mr = [r for r in min_rows if r["Case"] == case and r["Drug_ID"] not in EXCLUDED_IDS]
        ar = [r for r in auc_rows if r["Case"] == case and r["Drug_ID"] not in EXCLUDED_IDS]

        omin = [float(r["Original_Delta"]) for r in mr]
        vmin = [float(r["Validation_Delta_Active_minus_PBS"]) for r in mr]
        oauc = [float(r["Original_DeltaAUC"]) for r in ar]
        vauc = [float(r["Validation_DeltaAUC"]) for r in ar]

        rho_min, p_min = spearmanr(omin, vmin)
        rho_auc, p_auc = spearmanr(oauc, vauc)
        ci_min = bootstrap_spearman(omin, vmin, seed=101 + case_index)
        ci_auc = bootstrap_spearman(oauc, vauc, seed=201 + case_index)
        direction = np.mean([
            (a < 0 and b < 0) or (a > 0 and b > 0) or (a == 0 and b == 0)
            for a, b in zip(oauc, vauc)
        ])

        all_case = [r for r in min_rows if r["Case"] == case]
        active = np.asarray([float(r["Validation_Delta_Active_minus_PBS"]) for r in all_case])
        hi = np.asarray([float(r["Validation_Delta_HI_minus_PBS"]) for r in all_case])
        active_hi = active - hi
        test = wilcoxon(active, hi)
        specificity_p.append(test.pvalue)
        specificity.append([case, np.median(active_hi), test.pvalue, np.sum(active_hi <= -10)])

        summary.append([
            case, len(mr), rho_min, p_min, ci_min[0], ci_min[1],
            rho_auc, p_auc, ci_auc[0], ci_auc[1], direction
        ])

    qvalues = bh_fdr(specificity_p)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    with (out / "reproducibility_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "Case", "N", "DeltaMin_rho", "DeltaMin_p", "DeltaMin_CI_lower", "DeltaMin_CI_upper",
            "DeltaAUC_rho", "DeltaAUC_p", "DeltaAUC_CI_lower", "DeltaAUC_CI_upper",
            "DeltaAUC_direction_concordance"
        ])
        writer.writerows(summary)

    with (out / "specificity_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Case", "Median_Active_minus_HI", "Wilcoxon_p", "BH_FDR_q", "ActiveSpecific_n"])
        for row, q in zip(specificity, qvalues):
            writer.writerow(row[:3] + [q, row[3]])


if __name__ == "__main__":
    main()
