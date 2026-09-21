#!/usr/bin/env python3
"""
Compute the intrinsic pitch-angle profile of a GALFIT spiral model.

Supported GALFIT coordinate-rotation functions
-----------------------------------------------
    R0) power   -> alpha-tanh spiral
    R0) log     -> log-tanh spiral

For either model, the script reconstructs the GALFIT rotation law theta(r),
computes dtheta/dr analytically, and evaluates the local pitch angle

    p(r) = arctan[1 / |r dtheta/dr|],

where theta is converted to radians.

The returned pitch angle is intrinsic to the GALFIT spiral plane. R9
(inclination) and R10 (sky position angle) only project that intrinsic
spiral onto the sky, so no additional deprojection is applied here.

The script reports:
    - median pitch angle
    - 16th and 84th percentiles
    - asymmetric percentile uncertainties around the median
    - standard deviation
    - winding-reversal radii, if present

It also writes:
    - CSV radial profile
    - TXT statistical summary
    - PNG plot of pitch angle versus radius

Examples
--------
pitch_angle.py galfit.02

pitch_angle.py galfit.02 --component 2

pitch_angle.py galfit.02 \
    --median-rmin 30 --median-rmax 150

pitch_angle.py galfit.02 \
    --rmin 5 --rmax 220 \
    --median-rmin 30 --median-rmax 150 \
    --output-prefix ngc3627_pitch
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


CDEF = 0.23


def parse_galfit_file(filename: str):
    """Read GALFIT spiral parameters and plate scale."""
    text = Path(filename).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    components = {}
    current_component = None
    plate_scale = None

    component_re = re.compile(r"#\s*Component number:\s*(\d+)", re.IGNORECASE)
    rotation_re = re.compile(r"^\s*(R\d+)\)\s+([^\s#]+)", re.IGNORECASE)
    plate_re = re.compile(
        r"^\s*K\)\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
        r"\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
    )

    for line in lines:
        match = plate_re.match(line)
        if match:
            plate_scale = (float(match.group(1)), float(match.group(2)))

        match = component_re.search(line)
        if match:
            current_component = int(match.group(1))
            components.setdefault(current_component, {})
            continue

        if current_component is None:
            continue

        match = rotation_re.match(line)
        if not match:
            continue

        key = match.group(1).upper()
        raw_value = match.group(2)

        if key == "R0":
            components[current_component][key] = raw_value.lower()
        else:
            try:
                components[current_component][key] = float(raw_value)
            except ValueError:
                pass

    spiral_components = {
        number: pars
        for number, pars in components.items()
        if pars.get("R0", "none") in {"power", "log"}
    }

    return spiral_components, plate_scale


def tanh_transition(r, r_in, r_out, theta_out_deg):
    """
    GALFIT hyperbolic-tangent transition.

    The implementation follows the transition used by the GALFIT
    alpha-tanh and log-tanh coordinate-rotation models.
    """
    if r_out == r_in:
        raise ValueError("R1 and R2 cannot be equal.")
    if r_out == 0:
        raise ValueError("R2 cannot be zero.")

    a_value = 2.0 * CDEF / (abs(theta_out_deg) + CDEF) - 1.00001

    eps = 1.0e-12
    a_value = np.clip(a_value, -1.0 + eps, 1.0 - eps)

    b_value = (2.0 - np.arctanh(a_value)) * r_out / (r_out - r_in)

    u_value = b_value * (r / r_out - 1.0) + 2.0
    tanh_u = np.tanh(u_value)

    transition = 0.5 * (tanh_u + 1.0)

    # sech^2(u) = 1 - tanh^2(u), numerically safer than 1/cosh^2(u).
    dtransition_dr = 0.5 * (1.0 - tanh_u**2) * (b_value / r_out)

    return transition, dtransition_dr


def spiral_rotation_and_derivative(r, pars):
    """
    Return theta(r) in degrees and dtheta/dr in degrees per pixel.
    """
    mode = pars["R0"]
    r_in = pars["R1"]
    r_out = pars["R2"]
    theta_out = pars["R3"]
    shape_parameter = pars["R4"]

    transition, dtransition_dr = tanh_transition(
        r,
        r_in,
        r_out,
        theta_out,
    )

    theta = np.full_like(r, np.nan, dtype=float)
    dtheta_dr = np.full_like(r, np.nan, dtype=float)

    if mode == "power":
        alpha = shape_parameter

        base = 0.5 * (r / r_out + 1.0)
        valid = base > 0

        base_valid = base[valid]
        transition_valid = transition[valid]
        dtransition_valid = dtransition_dr[valid]

        power_term = base_valid**alpha
        dpower_dr = alpha * base_valid ** (alpha - 1.0) * (0.5 / r_out)

        theta[valid] = theta_out * transition_valid * power_term

        dtheta_dr[valid] = theta_out * (
            dtransition_valid * power_term + transition_valid * dpower_dr
        )

    elif mode == "log":
        r_ws = shape_parameter

        if r_ws == 0:
            raise ValueError("For R0=log, R4 (winding scale radius) cannot be zero.")

        normalization_argument = 1.0 + r_out / r_ws
        if normalization_argument <= 0:
            raise ValueError("Invalid log-tanh model: 1 + R2/R4 must be positive.")

        normalization = np.log(normalization_argument)
        if normalization == 0:
            raise ValueError("Invalid log-tanh normalization.")

        argument = 1.0 + r / r_ws
        valid = argument > 0

        radius_valid = r[valid]
        transition_valid = transition[valid]
        dtransition_valid = dtransition_dr[valid]

        log_term = np.log(1.0 + radius_valid / r_ws) / normalization

        dlog_dr = 1.0 / ((radius_valid + r_ws) * normalization)

        theta[valid] = theta_out * transition_valid * log_term

        dtheta_dr[valid] = theta_out * (
            dtransition_valid * log_term + transition_valid * dlog_dr
        )

    else:
        raise ValueError(f"Unsupported rotation function: {mode}")

    return theta, dtheta_dr


def pitch_angle(r, dtheta_dr_deg):
    """
    Compute the local intrinsic pitch angle in degrees.

    The astronomical pitch angle is the angle between the spiral tangent
    and the local circular direction.
    """
    dtheta_dr_rad = np.deg2rad(dtheta_dr_deg)
    winding_rate = np.abs(r * dtheta_dr_rad)

    # arctan2(1, 0) -> 90 deg, appropriate for a locally radial arm.
    pitch = np.rad2deg(
        np.arctan2(
            np.ones_like(winding_rate),
            winding_rate,
        )
    )

    pitch[~np.isfinite(dtheta_dr_deg)] = np.nan

    return pitch


def compute_statistics(r, pitch, radius_min, radius_max):
    """Compute pitch-angle statistics over a selected radial interval."""
    mask = (r >= radius_min) & (r <= radius_max) & np.isfinite(pitch)

    values = pitch[mask]

    if values.size == 0:
        raise ValueError(
            "No finite pitch-angle values in the requested statistics range."
        )

    p16, median, p84 = np.nanpercentile(
        values,
        [16.0, 50.0, 84.0],
    )

    std = np.nanstd(values, ddof=1) if values.size > 1 else 0.0

    return {
        "mask": mask,
        "n": int(values.size),
        "p16": float(p16),
        "median": float(median),
        "p84": float(p84),
        "lower_error": float(median - p16),
        "upper_error": float(p84 - median),
        "std": float(std),
    }


def find_reversal_radii(r, dtheta_dr):
    """Locate approximate radii where the winding direction reverses."""
    finite = np.isfinite(dtheta_dr)
    radius = r[finite]
    derivative = dtheta_dr[finite]

    if radius.size < 2:
        return []

    signs = np.sign(derivative)

    # Propagate the previous sign through isolated exact zeros.
    for index in range(1, signs.size):
        if signs[index] == 0:
            signs[index] = signs[index - 1]

    changes = np.where(signs[1:] * signs[:-1] < 0)[0]

    reversals = []

    for index in changes:
        r1 = radius[index]
        r2 = radius[index + 1]
        d1 = derivative[index]
        d2 = derivative[index + 1]

        if d2 != d1:
            root = r1 - d1 * (r2 - r1) / (d2 - d1)
        else:
            root = 0.5 * (r1 + r2)

        reversals.append(float(root))

    return reversals


def save_csv(
    filename,
    component,
    radius_pix,
    radius_arcsec,
    theta,
    dtheta_dr,
    pitch,
    stats_mask,
):
    """Write the full radial pitch-angle profile."""
    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8",
    ) as output:
        writer = csv.writer(output)

        writer.writerow(
            [
                "component",
                "radius_pix",
                "radius_arcsec",
                "theta_deg",
                "dtheta_dr_deg_per_pix",
                "pitch_deg",
                "in_statistics_range",
            ]
        )

        for values in zip(
            radius_pix,
            radius_arcsec,
            theta,
            dtheta_dr,
            pitch,
            stats_mask,
        ):
            writer.writerow(
                [
                    component,
                    values[0],
                    values[1],
                    values[2],
                    values[3],
                    values[4],
                    int(values[5]),
                ]
            )


def save_summary(
    filename,
    input_file,
    component,
    pars,
    stats,
    stats_rmin,
    stats_rmax,
    plate_scale,
    reversals,
):
    """Write a compact text summary."""
    lines = []

    lines.append(f"GALFIT file       : {input_file}")
    lines.append(f"Spiral component  : {component}")
    lines.append(f"Rotation function : {pars['R0']}")
    lines.append(f"R1 (inner radius) : {pars['R1']:.6g} pix")
    lines.append(f"R2 (outer radius) : {pars['R2']:.6g} pix")
    lines.append(f"R3 (rotation)     : {pars['R3']:.6g} deg")

    if pars["R0"] == "power":
        lines.append(f"R4 (alpha)        : {pars['R4']:.6g}")
    else:
        lines.append(f"R4 (r_ws)         : {pars['R4']:.6g} pix")

    if "R9" in pars:
        lines.append(f"R9 (inclination)  : {pars['R9']:.6g} deg")

    if "R10" in pars:
        lines.append(f"R10 (sky PA)      : {pars['R10']:.6g} deg")

    if plate_scale is not None:
        mean_scale = 0.5 * (plate_scale[0] + plate_scale[1])
        lines.append(f"Plate scale       : {mean_scale:.6g} arcsec/pix")

    lines.append(f"Statistics range  : {stats_rmin:.6g} -- {stats_rmax:.6g} pix")
    lines.append(f"Samples used      : {stats['n']}")
    lines.append(f"Pitch p16         : {stats['p16']:.6f} deg")
    lines.append(f"Pitch median      : {stats['median']:.6f} deg")
    lines.append(f"Pitch p84         : {stats['p84']:.6f} deg")
    lines.append(
        "Pitch median 16-84: "
        f"{stats['median']:.6f} "
        f"-{stats['lower_error']:.6f} "
        f"+{stats['upper_error']:.6f} deg"
    )
    lines.append(f"Pitch std         : {stats['std']:.6f} deg")

    if reversals:
        reversal_text = ", ".join(f"{value:.6g}" for value in reversals)
        lines.append(f"Winding reversal  : {reversal_text} pix")
    else:
        lines.append("Winding reversal  : none in sampled radial range")

    Path(filename).write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def save_plot(
    filename,
    plot_radius,
    pitch,
    stats,
    stats_xmin,
    stats_xmax,
    reversals_plot,
    radius_unit,
    component,
    mode,
):
    """Create pitch angle versus radius plot with the 16th-84th band."""
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for plotting. " "Install it or run with --no-plot."
        ) from exc

    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    ax.plot(
        plot_radius,
        pitch,
        linewidth=1.7,
        label="Local pitch angle",
    )

    # Shade the 16th-84th percentile interval only inside the
    # radial range used for the statistics.
    stats_x = (plot_radius >= stats_xmin) & (plot_radius <= stats_xmax)

    ax.fill_between(
        plot_radius[stats_x],
        stats["p16"],
        stats["p84"],
        alpha=0.20,
        label="16th-84th percentile",
    )

    ax.hlines(
        stats["median"],
        stats_xmin,
        stats_xmax,
        linestyles="--",
        linewidth=1.5,
        label=f"Median = {stats['median']:.2f} deg",
    )

    ax.axvline(
        stats_xmin,
        linestyle=":",
        linewidth=1.0,
    )
    ax.axvline(
        stats_xmax,
        linestyle=":",
        linewidth=1.0,
    )

    for reversal in reversals_plot:
        ax.axvline(
            reversal,
            linestyle="-.",
            linewidth=1.0,
            alpha=0.7,
        )

    summary_text = (
        f"$p_{{50}}$ = {stats['median']:.2f} deg\n"
        f"$p_{{16}}$ = {stats['p16']:.2f} deg\n"
        f"$p_{{84}}$ = {stats['p84']:.2f} deg\n"
        f"$\\sigma$ = {stats['std']:.2f} deg"
    )

    ax.text(
        0.97,
        0.97,
        summary_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.80,
        },
    )

    ax.set_xlabel(f"Radius [{radius_unit}]")
    ax.set_ylabel("Intrinsic pitch angle [deg]")
    ax.set_ylim(
        0.0,
        90.0,
    )

    ax.set_title(f"GALFIT pitch-angle profile: component {component} ({mode}-tanh)")

    ax.grid(
        alpha=0.25,
    )
    ax.legend(
        loc="best",
    )

    fig.tight_layout()
    fig.savefig(
        filename,
        dpi=200,
    )
    plt.close(fig)


def mainpitchAngle():
    parser = argparse.ArgumentParser(
        description=(
            "Compute intrinsic pitch angle versus radius from a "
            "GALFIT alpha-tanh or log-tanh spiral model."
        )
    )

    parser.add_argument(
        "galfit_file",
        help="GALFIT input/output parameter file.",
    )

    parser.add_argument(
        "-c",
        "--component",
        type=int,
        help=("Spiral component number. " "Default: first spiral component found."),
    )

    parser.add_argument(
        "--rmin",
        type=float,
        default=None,
        help=(
            "Minimum radius for the computed profile [pix]. "
            "Default: a small positive radius."
        ),
    )

    parser.add_argument(
        "--rmax",
        type=float,
        default=None,
        help=("Maximum radius for the computed profile [pix]. " "Default: R2."),
    )

    parser.add_argument(
        "--median-rmin",
        type=float,
        default=None,
        help=(
            "Minimum radius used for median and percentiles [pix]. "
            "Default: max(0, R1)."
        ),
    )

    parser.add_argument(
        "--median-rmax",
        type=float,
        default=None,
        help=("Maximum radius used for median and percentiles [pix]. " "Default: R2."),
    )

    parser.add_argument(
        "-n",
        "--npoints",
        type=int,
        default=2000,
        help="Number of radial samples. Default: 2000.",
    )

    parser.add_argument(
        "--output-prefix",
        default=None,
        help=("Output prefix. " "Default: <GALFIT filename>_pitch."),
    )

    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Do not create the PNG plot.",
    )

    args = parser.parse_args()

    spiral_components, plate_scale = parse_galfit_file(args.galfit_file)

    if not spiral_components:
        raise SystemExit(
            "No GALFIT spiral component with " "R0=power or R0=log was found."
        )

    if args.component is None:
        component = sorted(spiral_components)[0]
    else:
        component = args.component

        if component not in spiral_components:
            available = ", ".join(str(number) for number in sorted(spiral_components))
            raise SystemExit(
                f"Component {component} is not a spiral component. "
                f"Available spiral components: {available}"
            )

    pars = spiral_components[component]

    required = [
        "R0",
        "R1",
        "R2",
        "R3",
        "R4",
    ]

    missing = [key for key in required if key not in pars]

    if missing:
        raise SystemExit(f"Component {component} is missing: " f"{', '.join(missing)}")

    r_in = pars["R1"]
    r_out = pars["R2"]

    if r_out <= 0:
        raise SystemExit("R2 must be positive.")

    if args.npoints < 10:
        raise SystemExit("--npoints must be at least 10.")

    default_small_radius = max(
        1.0e-3,
        1.0e-4 * r_out,
    )

    radius_min = default_small_radius if args.rmin is None else args.rmin

    radius_max = r_out if args.rmax is None else args.rmax

    if radius_min < 0 or radius_max <= radius_min:
        raise SystemExit("Require 0 <= rmin < rmax.")

    radius = np.linspace(
        radius_min,
        radius_max,
        args.npoints,
    )

    theta, dtheta_dr = spiral_rotation_and_derivative(
        radius,
        pars,
    )

    pitch = pitch_angle(
        radius,
        dtheta_dr,
    )

    stats_rmin = (
        max(radius_min, max(0.0, r_in))
        if args.median_rmin is None
        else args.median_rmin
    )

    stats_rmax = (
        min(radius_max, r_out) if args.median_rmax is None else args.median_rmax
    )

    if stats_rmax <= stats_rmin:
        raise SystemExit("Require median-rmin < median-rmax.")

    stats = compute_statistics(
        radius,
        pitch,
        stats_rmin,
        stats_rmax,
    )

    reversals = find_reversal_radii(
        radius,
        dtheta_dr,
    )

    if plate_scale is not None:
        dx, dy = plate_scale
        mean_scale = 0.5 * (dx + dy)
        radius_arcsec = radius * mean_scale
    else:
        mean_scale = math.nan
        radius_arcsec = np.full_like(
            radius,
            np.nan,
        )

    input_path = Path(args.galfit_file)

    if args.output_prefix:
        prefix = Path(args.output_prefix)
    else:
        prefix = input_path.with_name(input_path.name + "_pitch")

    csv_file = str(prefix) + ".csv"
    txt_file = str(prefix) + "_summary.txt"
    png_file = str(prefix) + ".png"

    save_csv(
        csv_file,
        component,
        radius,
        radius_arcsec,
        theta,
        dtheta_dr,
        pitch,
        stats["mask"],
    )

    save_summary(
        txt_file,
        args.galfit_file,
        component,
        pars,
        stats,
        stats_rmin,
        stats_rmax,
        plate_scale,
        reversals,
    )

    if not args.no_plot:
        if np.isfinite(mean_scale):
            plot_radius = radius_arcsec
            stats_xmin = stats_rmin * mean_scale
            stats_xmax = stats_rmax * mean_scale
            reversals_plot = [value * mean_scale for value in reversals]
            radius_unit = "arcsec"
        else:
            plot_radius = radius
            stats_xmin = stats_rmin
            stats_xmax = stats_rmax
            reversals_plot = reversals
            radius_unit = "pix"

        save_plot(
            png_file,
            plot_radius,
            pitch,
            stats,
            stats_xmin,
            stats_xmax,
            reversals_plot,
            radius_unit,
            component,
            pars["R0"],
        )

    print(f"GALFIT file       : {args.galfit_file}")
    print(f"Spiral component  : {component}")
    print(f"Rotation function : {pars['R0']}")
    print(f"R1 (inner radius) : {pars['R1']:.6g} pix")
    print(f"R2 (outer radius) : {pars['R2']:.6g} pix")
    print(f"R3 (rotation)     : {pars['R3']:.6g} deg")

    if pars["R0"] == "power":
        print(f"R4 (alpha)        : {pars['R4']:.6g}")
    else:
        print(f"R4 (r_ws)         : {pars['R4']:.6g} pix")

    if "R9" in pars:
        print(f"R9 (inclination)  : {pars['R9']:.6g} deg")

    if "R10" in pars:
        print(f"R10 (sky PA)      : {pars['R10']:.6g} deg")

    if np.isfinite(mean_scale):
        print(f"Plate scale       : {mean_scale:.6g} arcsec/pix")

    print(f"Statistics range  : {stats_rmin:.6g} -- " f"{stats_rmax:.6g} pix")
    print(f"Samples used      : {stats['n']}")
    print(f"Pitch p16         : {stats['p16']:.6f} deg")
    print(f"Pitch median      : {stats['median']:.6f} deg")
    print(f"Pitch p84         : {stats['p84']:.6f} deg")
    print(
        "Pitch 16-84 range : "
        f"{stats['median']:.6f} "
        f"-{stats['lower_error']:.6f} "
        f"+{stats['upper_error']:.6f} deg"
    )
    print(f"Pitch std         : {stats['std']:.6f} deg")

    if reversals:
        formatted = ", ".join(f"{value:.3f}" for value in reversals)
        print("Winding reversal  : " f"{formatted} pix")
        print(
            "WARNING            : dtheta/dr changes sign. "
            "Interpret a single global pitch angle cautiously."
        )
    else:
        print("Winding reversal  : none in sampled range")

    print(f"CSV output        : {csv_file}")
    print(f"Summary output    : {txt_file}")

    if not args.no_plot:
        print(f"Plot output       : {png_file}")


if __name__ == "__main__":
    mainpitchAngle()
