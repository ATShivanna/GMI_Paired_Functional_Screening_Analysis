#!/usr/bin/env python3
"""Reverse-audit QC values from original raw 1536-well workbooks.

Raw clinical plate files are not included in the public release. This script documents the
audited calculation and can be run locally when authorized raw files are available.

High-response reference:
- PBS plate: Vehicle_DMSO plus PBS_only controls.
- GMI plate: Vehicle_DMSO controls.
Low-response reference:
- Lowest 5% of all valid raw plate values.

Zprime = 1 - 3 * (SD_high + SD_low) / (Mean_high - Mean_low)
Signal window = (Mean_high - Mean_low) / (SD_high + SD_low)
"""

from statistics import mean, stdev


def calculate_qc(high_values, all_plate_values):
    low_count = int(len(all_plate_values) * 0.05)
    low_values = sorted(all_plate_values)[:low_count]
    high_mean = mean(high_values)
    low_mean = mean(low_values)
    high_sd = stdev(high_values)
    low_sd = stdev(low_values)
    zprime = 1 - 3 * (high_sd + low_sd) / (high_mean - low_mean)
    signal_window = (high_mean - low_mean) / (high_sd + low_sd)
    return {
        "zprime": zprime,
        "signal_window": signal_window,
        "high_mean": high_mean,
        "low_mean": low_mean,
        "high_sd": high_sd,
        "low_sd": low_sd,
        "low_count": low_count,
    }
