#!/usr/bin/env python3

import argparse
import re
import os
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

    match = re.match(
        r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(.*)",
        value,
    )

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

    match = re.match(
        r"([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)(.*)",
        value,
    )

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


def get_transformed_parameters(match, axis_ratio=0.4, scale=1.0):
    """
    Compute the transformed ellipse parameters.

    Rules
    -----
    1. The center is unchanged.
    2. The first output ellipse is rotated by 90 degrees relative to the
       input major axis.
    3. The new major axis is the old minor axis.
    4. The new axis ratio is fixed by axis_ratio.
    5. The output axes are multiplied by scale.
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

    if axis_ratio <= 0 or axis_ratio > 1:
        raise ValueError("axis_ratio must satisfy 0 < axis_ratio <= 1.")

    if scale <= 0:
        raise ValueError("scale must be larger than 0.")

    # Identify the old major and minor semi-axes.
    if axis1 >= axis2:
        old_major = axis1
        old_minor = axis2
        old_major_angle = angle
    else:
        old_major = axis2
        old_minor = axis1
        old_major_angle = angle + 90.0

    # New major axis is the old minor axis.
    new_major = old_minor * scale

    # Fixed output axis ratio.
    new_minor = axis_ratio * new_major

    # First output ellipse: rotated 90 degrees from the input major axis.
    new_angle = (old_major_angle + 90.0) % 180.0

    return xcen, ycen, new_major, new_minor, new_angle, unit1


def transform_ellipse(match, axis_ratio=0.4, scale=1.0):
    """
    Transform one DS9 ellipse region.
    """
    xcen, ycen, new_major, new_minor, new_angle, unit = get_transformed_parameters(
        match,
        axis_ratio=axis_ratio,
        scale=scale,
    )

    return build_ellipse(
        xcen,
        ycen,
        new_major,
        new_minor,
        new_angle,
        unit,
    )


def transform_ellipse_rot90(match, axis_ratio=0.4, scale=1.0):
    """
    Transform one DS9 ellipse region and rotate it by an additional 90 degrees.
    """
    xcen, ycen, new_major, new_minor, new_angle, unit = get_transformed_parameters(
        match,
        axis_ratio=axis_ratio,
        scale=scale,
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


def transform_ellipse_rot45(match, axis_ratio=0.4, scale=1.0):
    """
    Transform one DS9 ellipse region and rotate it by an additional 45 degrees.
    """
    xcen, ycen, new_major, new_minor, new_angle, unit = get_transformed_parameters(
        match,
        axis_ratio=axis_ratio,
        scale=scale,
    )

    new_angle_rot45 = (new_angle + 45.0) % 180.0

    return build_ellipse(
        xcen,
        ycen,
        new_major,
        new_minor,
        new_angle_rot45,
        unit,
    )


def transform_region_file(
    input_file,
    output_file,
    output_file_rot90,
    output_file_rot45,
    axis_ratio=0.4,
    scale=1.0,
):
    """
    Read a DS9 region file and write two transformed region files.

    Parameters
    ----------
    input_file : str
        Input DS9 region file.
    output_file : str
        First output DS9 region file.
    output_file_rot90 : str
        Second output DS9 region file, rotated 90 degrees from the first.
    output_file_rot45 : str
        Third output DS9 region file, rotated 45 degrees from the first.
    axis_ratio : float, optional
        Axis ratio of the output ellipses. Default is 0.4.
    scale : float, optional
        Multiplicative scale factor for both output axes. Default is 1.0.
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    output_file_rot90 = Path(output_file_rot90)
    output_file_rot45 = Path(output_file_rot45)

    lines_out = []
    lines_out_rot90 = []
    lines_out_rot45 = []

    with input_file.open("r") as f:
        for line in f:
            if "ellipse" in line:
                new_line = ELLIPSE_RE.sub(
                    lambda match: transform_ellipse(
                        match,
                        axis_ratio=axis_ratio,
                        scale=scale,
                    ),
                    line,
                )

                new_line_rot90 = ELLIPSE_RE.sub(
                    lambda match: transform_ellipse_rot90(
                        match,
                        axis_ratio=axis_ratio,
                        scale=scale,
                    ),
                    line,
                )

                new_line_rot45 = ELLIPSE_RE.sub(
                    lambda match: transform_ellipse_rot45(
                        match,
                        axis_ratio=axis_ratio,
                        scale=scale,
                    ),
                    line,
                )

                lines_out.append(new_line)
                lines_out_rot90.append(new_line_rot90)
                lines_out_rot45.append(new_line_rot45)
            else:
                lines_out.append(line)
                lines_out_rot90.append(line)
                lines_out_rot45.append(line)

    with output_file.open("w") as f:
        f.writelines(lines_out)

    with output_file_rot90.open("w") as f:
        f.writelines(lines_out_rot90)

    with output_file_rot45.open("w") as f:
        f.writelines(lines_out_rot45)


def maintransformEllip():
    parser = argparse.ArgumentParser(
        description=(
            "Transform DS9 ellipse regions. The output ellipses use the "
            "old minor axis as the new major axis, use a fixed axis ratio, "
            "and can be enlarged or reduced with a scale factor."
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

    parser.add_argument(
        "--axis-ratio",
        type=float,
        default=0.4,
        help="Axis ratio b/a for the output ellipses. Default: 0.4.",
    )

    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help=(
            "Scale factor applied to both output axes. "
            "Use values larger than 1 to enlarge and smaller than 1 to shrink. "
            "Default: 1.0."
        ),
    )

    args = parser.parse_args()

    root, ext = os.path.splitext(args.output_region)

    output_region_rot90 = root + "_rot90" + ext
    output_region_rot45 = root + "_rot45" + ext

    transform_region_file(
        args.input_region,
        args.output_region,
        output_region_rot90,
        output_region_rot45,
        axis_ratio=args.axis_ratio,
        scale=args.scale,
    )


if __name__ == "__main__":
    maintransformEllip()
