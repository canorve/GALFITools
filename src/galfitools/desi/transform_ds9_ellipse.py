#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


ELLIPSE_RE = re.compile(
    r"ellipse\s*\(\s*"
    r"([^,]+)\s*,\s*"  # x center
    r"([^,]+)\s*,\s*"  # y center
    r"([^,]+)\s*,\s*"  # first axis
    r"([^,]+)\s*,\s*"  # second axis
    r"([^)]+)\s*"  # angle
    r"\)"
)


def parse_axis(value):
    """
    Parse a DS9 axis value.

    The function accepts values such as:
    20
    20.5
    20"
    3.2'

    Returns
    -------
    number : float
        Numerical value.
    unit : str
        Unit suffix, if present.
    """
    value = value.strip()

    match = re.match(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(.*)", value)

    if match is None:
        raise ValueError(f"Could not parse axis value: {value}")

    number = float(match.group(1))
    unit = match.group(2).strip()

    return number, unit


def parse_angle(value):
    """
    Parse a DS9 angle value in degrees.
    """
    value = value.strip()

    match = re.match(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(.*)", value)

    if match is None:
        raise ValueError(f"Could not parse angle value: {value}")

    angle = float(match.group(1))

    return angle


def format_value(value, unit=""):
    """
    Format numerical values for a DS9 region file.
    """
    return f"{value:.6f}{unit}"


def transform_ellipse(match):
    """
    Transform one DS9 ellipse region.

    Rules
    -----
    1. The center is unchanged.
    2. The position angle is rotated by 90 degrees.
    3. The new major axis is the old minor axis.
    4. The axis ratio is kept fixed.
    """
    xcen = match.group(1).strip()
    ycen = match.group(2).strip()

    axis1, unit1 = parse_axis(match.group(3))
    axis2, unit2 = parse_axis(match.group(4))
    angle = parse_angle(match.group(5))

    if unit1 != unit2:
        raise ValueError(
            "The two ellipse axes have different units: " f"{unit1!r} and {unit2!r}"
        )

    # Identify old major and minor semi-axes.
    if axis1 >= axis2:
        old_major = axis1
        old_minor = axis2
        old_major_angle = angle
    else:
        old_major = axis2
        old_minor = axis1
        old_major_angle = angle + 90.0

    old_axis_ratio = old_minor / old_major

    new_major = old_minor
    new_minor = new_major * old_axis_ratio

    new_angle = (old_major_angle + 90.0) % 180.0

    new_region = (
        f"ellipse("
        f"{xcen}, "
        f"{ycen}, "
        f"{format_value(new_major, unit1)}, "
        f"{format_value(new_minor, unit1)}, "
        f"{new_angle:.6f}"
        f")"
    )

    return new_region


def transform_region_file(input_file, output_file):
    """
    Read a DS9 region file and transform all ellipse regions.
    """
    input_file = Path(input_file)
    output_file = Path(output_file)

    new_lines = []

    with input_file.open("r") as f:
        for line in f:
            if "ellipse" in line:
                new_line = ELLIPSE_RE.sub(transform_ellipse, line)
                new_lines.append(new_line)
            else:
                new_lines.append(line)

    with output_file.open("w") as f:
        f.writelines(new_lines)


def maintransformEllip():
    parser = argparse.ArgumentParser(
        description=(
            "Transform DS9 ellipse regions. The new ellipse keeps the same "
            "center, rotates the position angle by 90 degrees, uses the old "
            "minor axis as the new major axis, and keeps the same axis ratio."
        )
    )

    parser.add_argument("input_region", help="Input DS9 region file.")

    parser.add_argument("output_region", help="Output DS9 region file.")

    args = parser.parse_args()

    transform_region_file(args.input_region, args.output_region)


if __name__ == "__main__":
    maintransformEllip()
