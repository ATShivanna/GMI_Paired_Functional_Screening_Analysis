import os
import pandas as pd
import numpy as np

BASE_PATH = "./data_raw/cases"

PRIMARY_CASES = ["A032", "A043", "A156", "A166", "A173", "A174", "A180"]
SECONDARY_CASES = ["A051", "A123", "A127"]

ALL_CASES = PRIMARY_CASES + SECONDARY_CASES

OUTDIR = "./normalized_output"
os.makedirs(OUTDIR, exist_ok=True)

def extract_plate_matrix(file_path):
    df = pd.read_excel(file_path, sheet_name="Plate results", header=None)

    start_row = None

    for i in range(len(df)):
        row = df.iloc[i, :10]

        # Count how many numeric values exist in first 10 columns
        numeric_count = pd.to_numeric(row, errors='coerce').notna().sum()

        # Detect header row (should contain many numeric values like 1–48)
        if numeric_count >= 5:
            start_row = i + 1
            break

    if start_row is None:
        raise ValueError(f"Plate start not found in {file_path}")

    # Extract 32 rows × 48 columns
    plate = df.iloc[start_row:start_row+32, 1:49]
    plate = plate.apply(pd.to_numeric, errors="coerce")

    if plate.shape != (32, 48):
        raise ValueError(f"Unexpected plate shape {plate.shape} in {file_path}")

    return plate

def get_vehicle_mask():
    # AB–AF rows, columns 41–48
    # row indices: A=0 ... Z=25, AA=26, AB=27, AC=28, AD=29, AE=30, AF=31
    rows = list(range(27, 32))   # AB–AF
    cols = list(range(40, 48))   # 41–48 in 0-based indexing

    mask = np.zeros((32, 48), dtype=bool)
    for r in rows:
        for c in cols:
            mask[r, c] = True
    return mask

def normalize_plate(plate, vehicle_mask):
    values = plate.values.astype(float)

    vehicle_vals = values[vehicle_mask]
    high = np.nanmedian(vehicle_vals)

    flat = values.flatten()
    flat = flat[~np.isnan(flat)]
    n_low = max(10, int(0.05 * len(flat)))
    low = np.nanmedian(np.sort(flat)[:n_low])

    denom = high - low
    if denom == 0:
        raise ValueError("Normalization denominator is zero")

    norm = (values - low) / denom
    norm = np.clip(norm, 0, 1) * 100

    return pd.DataFrame(norm)

vehicle_mask = get_vehicle_mask()
summary_rows = []

for case in ALL_CASES:
    case_group = "Primary" if case in PRIMARY_CASES else "Secondary"

    for condition in ["PBS", "GMI"]:
        file_path = f"{BASE_PATH}/{case}/T2_{condition}/{case}_T2_{condition}_3D.xlsx"

        if not os.path.exists(file_path):
            print(f"Missing file: {file_path}")
            summary_rows.append({
                "Case": case,
                "Condition": condition,
                "Group": case_group,
                "Status": "Missing file"
            })
            continue

        try:
            plate = extract_plate_matrix(file_path)
            norm_plate = normalize_plate(plate, vehicle_mask)

            out_file = os.path.join(
                OUTDIR,
                f"{case}_{condition}_normalized.csv"
            )
            norm_plate.to_csv(out_file, index=False, header=False)

            summary_rows.append({
                "Case": case,
                "Condition": condition,
                "Group": case_group,
                "Status": "OK",
                "Rows": norm_plate.shape[0],
                "Cols": norm_plate.shape[1],
                "Min": float(np.nanmin(norm_plate.values)),
                "Max": float(np.nanmax(norm_plate.values)),
                "Mean": float(np.nanmean(norm_plate.values))
            })

            print(f"{case} {condition} DONE | shape = {norm_plate.shape}")

        except Exception as e:
            print(f"{case} {condition} FAILED | {e}")
            summary_rows.append({
                "Case": case,
                "Condition": condition,
                "Group": case_group,
                "Status": f"FAILED: {e}"
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTDIR, "normalized_extraction_summary.csv"), index=False)

print("\nFinished.")
print(f"Output folder: {OUTDIR}")
print("Summary file: normalized_output/normalized_extraction_summary.csv")
