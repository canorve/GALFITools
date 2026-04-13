#!/usr/bin/env python3

"""
galfit_table_corrector.py

Usage example
-------------
    fitlogTableCorrector input.txt out.txt \
    --A 0.12 --K 0.04 \
    --pixscale 0.396 \
    --kpc-per-arcsec 0.48 \
    --re-units arcsec
"""
from __future__ import annotations
import argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Apply mag corrections and Re unit conversion to GALFIT fitlog table lines."
    )
    p.add_argument("input", type=Path, help="Input text file with GALFIT-style table.")
    p.add_argument("output", type=Path, help="Output text file.")
    # Photometric corrections
    p.add_argument(
        "--A",
        type=float,
        default=0.0,
        help="Extinction correction A_lambda (mag) to subtract.",
    )
    p.add_argument(
        "--K", type=float, default=0.0, help="K-correction (mag) to subtract."
    )
    # Geometry / scale
    p.add_argument(
        "--pixscale",
        type=float,
        default=None,
        help="Pixel scale in arcsec/pixel (required if --re-units is arcsec or kpc).",
    )
    p.add_argument(
        "--kpc-per-arcsec",
        type=float,
        default=None,
        help="Angular scale in kpc/arcsec (required if --re-units is kpc).",
    )
    p.add_argument(
        "--re-units",
        choices=["pix", "arcsec", "kpc"],
        default="arcsec",
        help="Unit for the effective radius (column 5) in the output.",
    )
    # Formatting options
    p.add_argument(
        "--mag-fmt", default=".2f", help="Format for magnitude column (e.g., .2f)."
    )
    p.add_argument("--re-fmt", default=".2f", help="Format for Re column (e.g., .2f).")
    p.add_argument(
        "--echo-original",
        action="store_true",
        help="Also emit the original component line as a comment just above the corrected one.",
    )
    return p.parse_args()


def ensure_scales(args):
    if args.re_units in ("arcsec", "kpc"):
        if args.pixscale is None:
            raise ValueError("--pixscale is required when --re-units is arcsec or kpc.")
    if args.re_units == "kpc":
        if args.kpc_per_arcsec is None:
            raise ValueError("--kpc-per-arcsec is required when --re-units is kpc.")


def convert_re(re_pix: float, args) -> float:
    if args.re_units == "pix":
        return re_pix
    re_arcsec = re_pix * args.pixscale
    if args.re_units == "arcsec":
        return re_arcsec
    # kpc
    return re_arcsec * args.kpc_per_arcsec


def format_number(x: float, spec: str) -> str:
    return f"{x:{spec}}"


def process_component_line(line: str, args) -> str:
    """
    Expected minimal structure:
    sersic, x, y, mag, re_pix, n, q, pa, <optional more numbers...>
    We correct col 4 (mag) and convert col 5 (Re) to requested units.
    All other columns are kept as-is.
    """
    # Split into comma-separated fields; strip whitespace per field.
    parts = [p.strip() for p in line.split(",")]
    if len(parts) < 7:
        # Not a proper component line; return unchanged.
        return line

    # Magnitude (col 4, index 3)
    try:
        mag = float(parts[3])
    except ValueError:
        return line  # leave untouched if unexpected

    mag_corr = mag - args.A - args.K
    parts[3] = format_number(mag_corr, args.mag_fmt)

    # Effective radius in pixels (col 5, index 4)
    try:
        re_pix = float(parts[4])
    except ValueError:
        return line

    re_new = convert_re(re_pix, args)
    parts[4] = format_number(re_new, args.re_fmt)

    # Reassemble with same comma-separated style
    return ",".join(parts)


def mainfitlogTableCorrector():
    args = parse_args()
    ensure_scales(args)

    txt_in = args.input.read_text(encoding="utf-8", errors="ignore").splitlines()
    out_lines = []
    for line in txt_in:
        # Identify model-component lines to modify. Here we target lines that start with 'sersic,'.
        if line.lstrip().startswith("sersic,"):
            if args.echo_original:
                out_lines.append("# ORIGINAL: " + line)
            fixed = process_component_line(line, args)
            out_lines.append(fixed)
        else:
            # Pass-through all non-component lines (headers, sky, chi^2, etc.)
            out_lines.append(line)

    args.output.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
