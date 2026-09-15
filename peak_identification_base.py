"""
Peak identification tool
Author: Rita Assis dos Santos (rita.assis.santos@tecnico.ulisboa.pt)
Copyright (c) 2026 Rita Assis dos Santos

Licensed under the MIT License. See the LICENSE file in the repository for the full license text.

Identify automatically GC peaks using a reference chromatogram.

Before running, update the CONFIGURATION section below.

Input: Tab-separated text files containing peak data (.txt files). Both the reference file and the files to be processed must be in the same folder as this script.

Returns: New tab-separated text files with the same name as the input files, but with "_identified" appended before the file extension.
The output files keep the columns "Peak #", "Sample", "Area S1", "Name", and "1st Dimension Time (min)" from the input files, but with the "Name" column updated to reflect the identified compounds.
"""

import glob
import os
from pathlib import Path
import pandas as pd

# ==================== CONFIGURATION ====================
# Edit this section only. The rest of the script can normally be left unchanged

# Input and output files
REFERENCE_FILE = "reference_chromatogram.txt"  # Tab-separated reference file reference file name (must be in the same folder as this script)
SCRIPT_DIR = Path(__file__).resolve().parent  # Folder containing the input files
OUTPUT_FOLDER = SCRIPT_DIR  # Folder where identified files are saved
FILE_PATTERN = (
    "*.txt"  # Input files to process (the reference is excluded): all .txt files
)
OUTPUT_SUFFIX = "_identified.txt"  # Added to each processed input filename
LIBRARY_PATTERN = "*.csv"  # Optional compound RT library in the script folder

# Column names in the input files
RT_COL = "1st Dimension Time (min)"  # Retention-time column
AREA_COL = "Area S1"  # Peak-area column
KEEP_COLUMNS = [
    "Peak #",
    "Sample",
    "Area S1",
    "Name",
    RT_COL,
]

# Peak identification
RT_TOLERANCE = 0.05  # Retention time tolerance in minutes for matching unknown peaks to reference peaks

# Post-processing thresholds
LOW_AREA_THRESHOLD = 20  # Remove every row with Area S1 below this value
EARLY_RT_THRESHOLD = 8  # Remove unidentified peaks at or below this RT

# Warning thresholds for unidentified peaks
WARNING_MEDIUM_MIN_AREA = 20
WARNING_MEDIUM_MAX_AREA = 100
WARNING_HIGH_AREA = 100

# Compounds whose rows should be removed during post-processing
COMPOUNDS_TO_REMOVE = [
    "bis(trimethylsilyl)trifluoroacetamide",
    "tetrasiloxane, decamethyl-",
    "column bleeding",
    "methylamine",
    "trimethylsilyl",
    "_der",
    "di-trimethylsilyl peroxide",
]


# ==================== HELPERS ====================


def load_raw(path):
    """
    Load a raw tab-separated GC-GC file and keep only the needed columns.
    Names of the columns required: "Peak #", "Sample", "Area S1", "Name", "1st Dimension Time (min)".
    Convert "Area S1" and "1st Dimension Time (min)" to numeric.

    Returns a dataframe with the same columns, or raises ValueError if any are missing.
    """
    df = pd.read_csv(path, sep="\t")

    missing = [column for column in KEEP_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing expected column(s): {missing}")

    df = df[KEEP_COLUMNS].copy()
    df[AREA_COL] = pd.to_numeric(df[AREA_COL], errors="coerce")
    df[RT_COL] = pd.to_numeric(df[RT_COL], errors="coerce")
    return df


def clean_df(df):
    """
    - Drop rows where "Peak #" is empty.
    - Remove duplicates (same "Peak #").

    Returns a cleaned dataframe with the same columns.
    """
    df = df.copy()
    df = df[df["Peak #"].notna()]
    return df.drop_duplicates(subset=["Peak #"], keep="first")


def is_unknown_peak(name):
    """
    True for unlabeled rows such as "peak", "peak 1", "peak 1:2".
    """
    if not isinstance(name, str):
        return False
    return "peak" in name.lower()


# ==================== REFERENCE BUILDING ====================


def print_unknown_warnings(df, heading):
    """
    Print configured warnings for unidentified peaks in a dataframe.
    """
    unknown_mask = df["Name"].apply(is_unknown_peak)
    high_mask = unknown_mask & (df[AREA_COL] > WARNING_HIGH_AREA)
    medium_mask = unknown_mask & df[AREA_COL].between(
        WARNING_MEDIUM_MIN_AREA, WARNING_MEDIUM_MAX_AREA
    )

    if high_mask.any():
        print(
            f"\nMASTER WARNING{heading}: Unidentified peaks with {AREA_COL} > "
            f"{WARNING_HIGH_AREA}:"
        )
        print(df.loc[high_mask, ["Peak #", "Name", AREA_COL, RT_COL]])

    if medium_mask.any():
        print(
            f"\nWarning{heading}: Unidentified peaks with "
            f"{WARNING_MEDIUM_MIN_AREA} <= {AREA_COL} <= "
            f"{WARNING_MEDIUM_MAX_AREA}:"
        )
        print(df.loc[medium_mask, ["Peak #", "Name", AREA_COL, RT_COL]])


def build_reference(df_ref):
    """
    Build a reference list from identified, valid reference peaks.

    Returns a list of dictionaries with keys "compound" and "rt".
    """
    df = clean_df(df_ref)
    print_unknown_warnings(df, "")

    before = len(df)
    df = df[~df["Name"].apply(is_unknown_peak)]
    print(f"Reference: removed {before - len(df)} unidentified row(s)")

    df = df[df[RT_COL].notna()]
    reference = [
        {"compound": row["Name"], "rt": float(row[RT_COL])} for _, row in df.iterrows()
    ]

    print(f"Reference library built: {len(reference)} identified peak(s)\n")
    return reference


# ==================== IDENTIFICATION ====================


def identify_peaks(df, reference, rt_col=RT_COL, tolerance=RT_TOLERANCE):
    """
    For every unidentified row in df, look for reference compounds
    whose RT window contains the peak's RT. Each reference peak is considered
    independently, including duplicate compound names. If several match, keep
    the closest.

    Returns a results dataframe (one row per unidentified peak).
    """
    results = []

    for index, row in df.iterrows():
        name = row["Name"]
        if not is_unknown_peak(name):
            continue

        peak_rt = row[rt_col]
        if pd.isna(peak_rt):
            continue

        matches = []
        for reference_peak in reference:
            distance = abs(peak_rt - reference_peak["rt"])
            if distance <= tolerance:
                matches.append(
                    {
                        "compound": reference_peak["compound"],
                        "distance": distance,
                    }
                )

        if matches:
            best = min(matches, key=lambda match: match["distance"])
            identified_as = best["compound"]
            all_matches = [match["compound"] for match in matches]
        else:
            identified_as = name
            all_matches = []

        results.append(
            {
                "idx": index,
                "original_name": name,
                "retention_time": peak_rt,
                "identified_as": identified_as,
                "all_matches": all_matches,
            }
        )

    return pd.DataFrame(results)


def apply_identifications(df, identification_results):
    """
    Update the "Name" column of df with the identified compound names from identification_results.

    Returns a new dataframe with the same columns as df.
    """
    df_out = df.copy()
    for _, result in identification_results.iterrows():
        if result["idx"] in df_out.index:
            df_out.at[result["idx"], "Name"] = result["identified_as"]
    return df_out


def load_rt_library(path):
    """
    Load an optional compound RT library and calculate its matching statistics.

    The library must contain "Name", "Min_RT_min", and "Max_RT_min". The mean
    is the midpoint of the interval and the standard deviation is the
    population standard deviation of the two interval endpoints, which is
    equivalent to half the interval width.
    """
    df = pd.read_csv(path)
    required_columns = ["Name", "Min_RT_min", "Max_RT_min"]
    missing = [column for column in required_columns if column not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing expected column(s): {missing}")

    df = df[required_columns].copy()
    df["Name"] = df["Name"].astype(str).str.strip()
    df["Min_RT_min"] = pd.to_numeric(df["Min_RT_min"], errors="coerce")
    df["Max_RT_min"] = pd.to_numeric(df["Max_RT_min"], errors="coerce")
    df = df.dropna(subset=["Name", "Min_RT_min", "Max_RT_min"])
    df = df[df["Name"] != ""]

    lower = df[["Min_RT_min", "Max_RT_min"]].min(axis=1)
    upper = df[["Min_RT_min", "Max_RT_min"]].max(axis=1)
    df["Min_RT_min"] = lower
    df["Max_RT_min"] = upper
    df["Mean_RT_min"] = (lower + upper) / 2
    df["Std_RT_min"] = (upper - lower) / 2
    return df


def identify_peaks_from_library(df, library, rt_col=RT_COL):
    """
    Identify only still-unknown peaks using compound RT intervals.

    If intervals overlap, the compound with the closest mean RT is selected.
    """
    results = []
    for index, row in df.iterrows():
        if not is_unknown_peak(row["Name"]) or pd.isna(row[rt_col]):
            continue

        matches = []
        for _, library_row in library.iterrows():
            match_min = library_row["Mean_RT_min"] - library_row["Std_RT_min"]
            match_max = library_row["Mean_RT_min"] + library_row["Std_RT_min"]
            if match_min <= row[rt_col] <= match_max:
                matches.append(
                    {
                        "compound": library_row["Name"],
                        "distance": abs(row[rt_col] - library_row["Mean_RT_min"]),
                    }
                )

        if matches:
            best = min(matches, key=lambda match: match["distance"])
            identified_as = best["compound"]
            all_matches = [match["compound"] for match in matches]
        else:
            identified_as = row["Name"]
            all_matches = []

        results.append(
            {
                "idx": index,
                "original_name": row["Name"],
                "retention_time": row[rt_col],
                "identified_as": identified_as,
                "all_matches": all_matches,
            }
        )

    return pd.DataFrame(results)


def find_rt_library():
    """Return the first valid CSV RT library in the script folder, if any."""
    required_columns = {"Name", "Min_RT_min", "Max_RT_min"}
    for path in sorted(SCRIPT_DIR.glob(LIBRARY_PATTERN)):
        try:
            columns = set(pd.read_csv(path, nrows=0).columns)
        except (OSError, pd.errors.ParserError):
            continue
        if required_columns.issubset(columns):
            return path
    return None


# ==================== POST-PROCESSING ====================


def compute_unknown_percentage(df):
    """
    Percentage of total Area S1 that belongs to unidentified ('peak'/'ni') rows.

    Returns a float between 0 and 100.
    """
    total_area = df[AREA_COL].sum()
    if not total_area:
        return 0.0

    peak_mask = df["Name"].apply(is_unknown_peak)
    names = df["Name"].astype(str).str.strip().str.lower()
    ni_mask = names.str.match(r"^ni(?![a-z])")
    unknown_area = df.loc[peak_mask | ni_mask, AREA_COL].sum()
    return unknown_area / total_area * 100


def post_process_df(df, rt_col=RT_COL):
    """
    Apply filtering without cleaning or normalizing compound names.

    Returns a tuple (df_cleaned, unknown_percentage).
    """
    df = df.copy()

    print(f"\nRemoving low area peaks ({AREA_COL} < {LOW_AREA_THRESHOLD})...")
    low_area_mask = df[AREA_COL] < LOW_AREA_THRESHOLD
    print(f"Number of peaks removed: {low_area_mask.sum()}")
    df = df[~low_area_mask]
    print(f"Remaining rows after low area removal: {len(df)}")

    print(f"\nRemoving unidentified peaks with R.T. <= " f"{EARLY_RT_THRESHOLD} min...")
    early_peak_mask = df["Name"].str.lower().str.contains("peak", na=False) & (
        df[rt_col] <= EARLY_RT_THRESHOLD
    )
    print(f"Number of early 'peak' rows removed: {early_peak_mask.sum()}")
    df = df[~early_peak_mask]
    print(f"Remaining rows after early peak removal: {len(df)}")

    print("\nRemoving specific unwanted compounds...")
    mask = pd.Series(False, index=df.index)
    for compound in COMPOUNDS_TO_REMOVE:
        mask |= (
            df["Name"].str.lower().str.contains(compound.lower(), na=False, regex=False)
        )

    removed = df[mask]
    print(f"Compounds being removed: {len(removed)}")
    if len(removed) > 0:
        print("Removed compound names:")
        print(removed[["Name", rt_col, AREA_COL]].to_string())

    df = df[~mask]
    print(f"Remaining rows after compound removal: {len(df)}")

    print_unknown_warnings(df, " (post-processing)")
    unknown_pct = compute_unknown_percentage(df)
    print(f"\nUnknown peak area: {unknown_pct:.2f}% of total area")
    return df, unknown_pct


# ==================== PIPELINE ====================


def process_file(path, reference, rt_library=None):
    print("=" * 100)
    print(f"Processing: {path}")
    df = clean_df(load_raw(path))
    print(f"Loaded {len(df)} rows")

    identification_results = identify_peaks(df, reference)
    n_identified = (
        (
            identification_results["identified_as"]
            != identification_results["original_name"]
        ).sum()
        if len(identification_results)
        else 0
    )
    n_ambiguous = (
        (identification_results["all_matches"].apply(len) > 1).sum()
        if len(identification_results)
        else 0
    )

    df_out = apply_identifications(df, identification_results)
    print(f"Unidentified rows found: {len(identification_results)}")
    print(f"  -> newly identified: {n_identified}")
    print(f"  -> still unknown: {len(identification_results) - n_identified}")
    print(f"  -> ambiguous matches: {n_ambiguous}")

    if rt_library is not None:
        library_results = identify_peaks_from_library(df_out, rt_library)
        library_identified = (
            (library_results["identified_as"] != library_results["original_name"]).sum()
            if len(library_results)
            else 0
        )
        df_out = apply_identifications(df_out, library_results)
        print(
            f"CSV library pass: identified {library_identified} of "
            f"{len(library_results)} still-unknown row(s)"
        )

    return df_out


def main():
    print("Loading reference file...")
    reference_path = (SCRIPT_DIR / REFERENCE_FILE).resolve()
    reference = build_reference(load_raw(reference_path))

    library_path = find_rt_library()
    rt_library = None
    if library_path is None:
        print("No CSV RT library found. Skipping the optional library pass.")
    else:
        try:
            rt_library = load_rt_library(library_path)
            print(
                f"Loaded CSV RT library: {library_path} ({len(rt_library)} compound range(s))"
            )
        except (OSError, ValueError, pd.errors.ParserError) as error:
            print(f"Warning: Could not load CSV RT library: {error}")

    files = sorted(
        path
        for path in glob.glob(os.path.join(SCRIPT_DIR, FILE_PATTERN))
        if Path(path).resolve() != reference_path
        and not path.endswith("_identified.txt")
    )

    if not files:
        print(f"No input files matching {FILE_PATTERN} found in '{SCRIPT_DIR}'.")
        return

    summary = []
    for file_path in files:
        df_out = process_file(file_path, reference, rt_library)
        df_clean, unknown_pct = post_process_df(df_out)

        output_name = os.path.splitext(os.path.basename(file_path))[0] + OUTPUT_SUFFIX
        output_path = OUTPUT_FOLDER / output_name
        df_clean.to_csv(output_path, sep="\t", index=False)
        print(f"Saved -> {output_path}")
        summary.append((os.path.basename(file_path), unknown_pct))

    print("\n" + "=" * 100)
    print(f"DONE. Processed {len(files)} file(s). Outputs in '{OUTPUT_FOLDER}/'.")
    print("\nUnknown peak area (% of total) by file:")
    for name, percentage in summary:
        print(f"  {name}: {percentage:.2f}%")
    print("=" * 100)


if __name__ == "__main__":
    main()
