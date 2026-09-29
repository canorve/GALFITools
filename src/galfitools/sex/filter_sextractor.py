#!/usr/bin/env python3

"""Filter 15-column SExtractor ASCII or ASCII_HEAD catalogs.

Usage:
    filterSex hot.cat filtered.cat
    filterSex hot.cat filtered.cat --flags 4 --mag 18 --class-star 0.6

All three upper limits are exclusive. No third-party packages are required.
ASCII_HEAD columns are located by their header names. Headerless ASCII must
use the example's order:
NUMBER ALPHA_J2000 DELTA_J2000 XPEAK_IMAGE YPEAK_IMAGE MAG_BEST KRON_RADIUS
FLUX_RADIUS ISOAREA_IMAGE A_IMAGE ELLIPTICITY THETA_IMAGE BACKGROUND CLASS_STAR FLAGS

Retained lines and comments are copied verbatim; objects are not renumbered.
Nonfinite filter values are rejected. Invalid row widths or nonnumeric data
raise an error before an output file is opened.
"""

import math
from pathlib import Path
import re


COLUMN_COUNT = 15
DEFAULT_COLUMNS = {"MAG_BEST": 5, "CLASS_STAR": 13, "FLAGS": 14}
HEADER_PATTERN = re.compile(r"^\s*#\s*(\d+)\s+(\S+)")


def filter_catalog(
    input_path,
    output_path,
    flags_max=4,
    mag_max=18.0,
    class_star_max=0.6,
    overwrite=False,
):
    """Write matching rows; return (input_count, retained_count).

    ASCII_HEAD must declare each of the 15 scalar columns exactly once.
    The output format is inherited from the input.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output must be different files.")
    if output_path.exists() and input_path.samefile(output_path):
        raise ValueError("Input and output refer to the same file.")
    if not all(math.isfinite(value) for value in (flags_max, mag_max, class_star_max)):
        raise ValueError("All thresholds must be finite.")

    with input_path.open(encoding="utf-8", newline="") as source:
        lines = source.readlines()

    headers = {}
    for line_number, line in enumerate(lines, 1):
        match = HEADER_PATTERN.match(line)
        if match:
            column = int(match.group(1)) - 1
            if column in headers:
                raise ValueError(f"Line {line_number}: duplicate column definition.")
            headers[column] = match.group(2)

    columns = DEFAULT_COLUMNS.copy()
    if headers:
        if set(headers) != set(range(COLUMN_COUNT)):
            raise ValueError("ASCII_HEAD must define exactly columns 1 through 15.")
        for name in columns:
            matches = [index for index, label in headers.items() if label == name]
            if len(matches) != 1:
                raise ValueError(f"Header must contain exactly one {name} column.")
            columns[name] = matches[0]

    output_lines = []
    total = kept = 0
    for line_number, line in enumerate(lines, 1):
        if not line.strip() or line.lstrip().startswith("#"):
            output_lines.append(line)
            continue
        # Allow an optional trailing comment without counting it as a column.
        fields = line.split("#", 1)[0].split()
        if len(fields) != COLUMN_COUNT:
            raise ValueError(
                f"Line {line_number}: expected 15 columns, found {len(fields)}."
            )
        try:
            values = [float(field) for field in fields]
        except ValueError as error:
            raise ValueError(
                f"Line {line_number}: nonnumeric catalog value."
            ) from error
        flags = values[columns["FLAGS"]]
        magnitude = values[columns["MAG_BEST"]]
        class_star = values[columns["CLASS_STAR"]]
        total += 1
        if not all(math.isfinite(value) for value in (flags, magnitude, class_star)):
            continue
        if flags < 0 or not flags.is_integer():
            raise ValueError(
                f"Line {line_number}: FLAGS must be a nonnegative integer."
            )
        if flags < flags_max and magnitude < mag_max and class_star < class_star_max:
            output_lines.append(line)
            kept += 1

    # Exclusive creation protects existing files unless overwrite is requested.
    with output_path.open(
        "w" if overwrite else "x", encoding="utf-8", newline=""
    ) as target:
        target.writelines(output_lines)
    return total, kept
