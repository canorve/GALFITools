#!/usr/bin/env python3
"""
Read a list of GALFIT parameter files, extract model parameters,
compute statistics across files, and save the results to a CSV file.

Assumptions
-----------
- All GALFIT files correspond to the same model.
- Components appear in the same order in all files.
- The uncertainty of each parameter is estimated from the scatter
  among the different GALFIT solutions.

Usage
-----
galfit_stats_to_csv.py galfit_list.txt output.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from statistics import mean, median, stdev
from typing import Dict, List, Tuple, Any


PARAMETER_LABELS = {
    "sersic": {
        "1": "x_center,y_center",
        "3": "magnitude",
        "4": "r_e",
        "5": "sersic_n",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "expdisk": {
        "1": "x_center,y_center",
        "3": "magnitude",
        "4": "scale_length",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "devauc": {
        "1": "x_center,y_center",
        "3": "magnitude",
        "4": "r_e",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "gaussian": {
        "1": "x_center,y_center",
        "3": "magnitude",
        "4": "fwhm",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "ferrer": {
        "1": "x_center,y_center",
        "3": "central_surface_brightness_or_mag",
        "4": "outer_truncation_radius",
        "5": "alpha",
        "6": "beta",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "moffat": {
        "1": "x_center,y_center",
        "3": "magnitude",
        "4": "fwhm",
        "5": "powerlaw",
        "9": "axis_ratio",
        "10": "position_angle",
    },
    "sky": {
        "1": "sky_background",
        "2": "dsky_dx",
        "3": "dsky_dy",
    },
}


def is_float(token: str) -> bool:
    """Return True if token can be interpreted as a float."""
    try:
        float(token)
        return True
    except ValueError:
        return False


def clean_line(line: str) -> str:
    """Remove comments and surrounding blanks."""
    if "#" in line:
        line = line.split("#", 1)[0]
    return line.strip()


def extract_numeric_values(rest: str) -> List[float]:
    """
    Extract numeric values from the right-hand side of a GALFIT line.

    This ignores fit flags such as trailing '1 1' or '0'.
    It only keeps tokens that can be converted to float.
    """
    tokens = rest.replace(",", " ").split()
    values = [float(token) for token in tokens if is_float(token)]
    return values


def parse_galfit_file(filepath: Path) -> List[Dict[str, Any]]:
    """
    Parse one GALFIT parameter file.

    Returns
    -------
    components : list of dict
        Each element contains one parameter entry:
        {
            "component_index": int,
            "component_type": str,
            "parameter_number": str,
            "parameter_label": str,
            "value_index": int,
            "value": float
        }
    """
    results: List[Dict[str, Any]] = []
    component_index = -1
    component_type = "unknown"

    with filepath.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = clean_line(raw_line)
            if not line:
                continue

            match = re.match(r"^([A-Za-z0-9]+)\)\s*(.*)$", line)
            if not match:
                continue

            parameter_number = match.group(1)
            rest = match.group(2).strip()

            if parameter_number == "0":
                component_index += 1
                component_type = rest.split()[0].lower() if rest else "unknown"
                continue

            if component_index < 0:
                continue
            if parameter_number.upper() == "Z":
                continue

            values = extract_numeric_values(rest)
            if not values:
                continue

            parameter_label = PARAMETER_LABELS.get(component_type, {}).get(
                parameter_number, f"param_{parameter_number}"
            )

            for value_index, value in enumerate(values):
                sub_label = parameter_label
                if len(values) == 2 and parameter_label == "x_center,y_center":
                    sub_label = "x_center" if value_index == 0 else "y_center"
                elif len(values) > 1:
                    sub_label = f"{parameter_label}_{value_index + 1}"

                results.append(
                    {
                        "component_index": component_index,
                        "component_type": component_type,
                        "parameter_number": parameter_number,
                        "parameter_label": sub_label,
                        "value_index": value_index,
                        "value": value,
                    }
                )

    return results


def build_statistics(
    all_data: Dict[Tuple[int, str, str, str, int], List[float]]
) -> List[Dict[str, Any]]:
    """Compute statistics for each grouped parameter."""
    rows: List[Dict[str, Any]] = []

    for key, values in sorted(all_data.items()):
        (
            component_index,
            component_type,
            parameter_number,
            parameter_label,
            value_index,
        ) = key

        n_values = len(values)

        if n_values == 1:
            std_value = 0.0
            stderr_value = 0.0
        else:
            std_value = stdev(values)
            stderr_value = std_value / math.sqrt(n_values)

        rows.append(
            {
                "component_index": component_index,
                "component_type": component_type,
                "parameter_number": parameter_number,
                "parameter_label": parameter_label,
                "value_index": value_index,
                "n_files": n_values,
                "mean": mean(values),
                "stddev": std_value,
                "stderr_mean": stderr_value,
                "median": median(values),
                "min": min(values),
                "max": max(values),
            }
        )

    return rows


def read_file_list(list_file: Path) -> List[Path]:
    """Read a text file containing GALFIT file paths, one per line."""
    paths: List[Path] = []

    with list_file.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            paths.append(Path(line))

    return paths


def write_csv(rows: List[Dict[str, Any]], output_csv: Path) -> None:
    """Write statistics rows to CSV."""
    fieldnames = [
        "component_index",
        "component_type",
        "parameter_number",
        "parameter_label",
        "value_index",
        "n_files",
        "mean",
        "stddev",
        "stderr_mean",
        "median",
        "min",
        "max",
    ]

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def mainGalStats() -> None:
    parser = argparse.ArgumentParser(
        description="Compute GALFIT parameter statistics from multiple parameter files."
    )
    parser.add_argument(
        "list_file",
        type=Path,
        help="Text file with the list of GALFIT parameter files, one per line.",
    )
    parser.add_argument(
        "output_csv",
        type=Path,
        help="Output CSV filename.",
    )
    args = parser.parse_args()

    galfit_files = read_file_list(args.list_file)

    if not galfit_files:
        raise SystemExit("Error: the list file is empty or contains no valid entries.")

    grouped_data: Dict[Tuple[int, str, str, str, int], List[float]] = {}

    for galfit_file in galfit_files:
        if not galfit_file.exists():
            print(f"Warning: file not found, skipped: {galfit_file}")
            continue

        parsed = parse_galfit_file(galfit_file)

        for item in parsed:
            key = (
                item["component_index"],
                item["component_type"],
                item["parameter_number"],
                item["parameter_label"],
                item["value_index"],
            )
            grouped_data.setdefault(key, []).append(item["value"])

    if not grouped_data:
        raise SystemExit("Error: no GALFIT parameters could be parsed.")

    rows = build_statistics(grouped_data)
    write_csv(rows, args.output_csv)

    print(f"Done. CSV written to: {args.output_csv}")
