#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
import os

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

    Parameters
    ----------
    value : str
        Axis value from a DS9 ellipse region.

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

    Parameters
    ----------
    value : str
        Angle value from a DS9 ellipse region.

    Returns
    -------
    angle : float
        Angle in degrees.
    """
    value = value.strip()

    match = re.match(r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(.*)", value)

    if match is None:
        raise ValueError(f"Could not parse angle value: {value}")

    angle = float(match.group(1))

    return angle


def format_value(value, unit=""):
    """
    Format a numerical value for a DS9 region file.
    """
    return f"{value:.6f}{unit}"


def build_ellipse(xcen, ycen, major_axis, minor_axis, angle, unit):
    """
    Build a DS9 ellipse string.
    """
    return (
        f"ellipse("
        f"{xcen}, "
        f"{ycen}, "
        f"{format_value(major_axis, unit)}, "
        f"{format_value(minor_axis, unit)}, "
        f"{angle:.6f}"
        f")"
    )


def get_transformed_parameters(match):
    """
    Compute the transformed ellipse parameters.

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

    return xcen, ycen, new_major, new_minor, new_angle, unit1


def transform_ellipse(match):
    """
    Transform one DS9 ellipse region.
    """
    xcen, ycen, new_major, new_minor, new_angle, unit = get_transformed_parameters(
        match
    )

    return build_ellipse(
        xcen,
        ycen,
        new_major,
        new_minor,
        new_angle,
        unit,
    )


def transform_ellipse_rot90(match):
    """
    Transform one DS9 ellipse region and rotate it by an additional 90 degrees.
    """
    xcen, ycen, new_major, new_minor, new_angle, unit = get_transformed_parameters(
        match
    )

    new_angle_rot90 = (new_angle + 90.0) % 180.0

    return build_ellipse(
        xcen,
        ycen,
        new_major,
        new_minor,
        new_angle_rot90,
        unit,
    )


def transform_region_file(input_file, output_file, output_file_rot90):
    """
    Read a DS9 region file and write two transformed region files.

    Parameters
    ----------
    input_file : str
        Input DS9 region file.
    output_file : str
        Output DS9 region file with the transformed ellipse.
    output_file_rot90 : str
        Output DS9 region file with the transformed ellipse rotated by 90 degrees.
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    output_file_rot90 = Path(output_file_rot90)

    lines_out = []
    lines_out_rot90 = []

    with input_file.open("r") as f:
        for line in f:
            if "ellipse" in line:
                new_line = ELLIPSE_RE.sub(transform_ellipse, line)
                new_line_rot90 = ELLIPSE_RE.sub(transform_ellipse_rot90, line)

                lines_out.append(new_line)
                lines_out_rot90.append(new_line_rot90)
            else:
                lines_out.append(line)
                lines_out_rot90.append(line)

    with output_file.open("w") as f:
        f.writelines(lines_out)

    with output_file_rot90.open("w") as f:
        f.writelines(lines_out_rot90)


def maintransformEllip():
    parser = argparse.ArgumentParser(
        description=(
            "Transform DS9 ellipse regions. The new ellipse keeps the same "
            "center, rotates the position angle by 90 degrees, uses the old "
            "minor axis as the new major axis, and keeps the same axis ratio."
            "It also creates a second output with the same axes "
            "but rotated by an additional 90 degrees."
        )
    )

    parser.add_argument(
        "input_region",
        help="Input DS9 region file.",
    )

    parser.add_argument(
        "output_region",
        help="First output DS9 region file.",
    )

    args = parser.parse_args()

    root, extension = os.path.splitext(args.output_region)

    output_region_rot90 = root + "_rot90" + extension

    transform_region_file(
        args.input_region,
        args.output_region,
        output_region_rot90,
    )


if __name__ == "__main__":
    maintransformEllip()
