#!/usr/bin/env python3
"""Generate randomized GALFIT initial files by calling getSersic."""

from __future__ import annotations

import argparse
import random
import re
import tempfile
from pathlib import Path
from typing import Callable


FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
ELLIPSE_PATTERN = re.compile(
    rf"\bellipse\s*\(\s*"
    rf"(?P<x>{FLOAT_PATTERN})\s*,\s*"
    rf"(?P<y>{FLOAT_PATTERN})\s*,\s*"
    rf"(?P<rx>{FLOAT_PATTERN})\s*,\s*"
    rf"(?P<ry>{FLOAT_PATTERN})\s*,\s*"
    rf"(?P<angle>{FLOAT_PATTERN})\s*\)",
    re.IGNORECASE,
)
COORDINATE_PATTERN = re.compile(
    r"\b(physical|image|fk4|fk5|icrs|galactic|ecliptic|linear|amplifier|detector)\b",
    re.IGNORECASE,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate randomized GALFIT parameter files by repeatedly calling "
            "the getSersic function."
        )
    )
    parser.add_argument("Image", help="input FITS image")
    parser.add_argument("RegFile", help="DS9 ellipse region for the galaxy")
    parser.add_argument(
        "-N",
        "--number",
        type=positive_integer,
        required=True,
        help="number of GALFIT initial-parameter files to generate",
    )
    parser.add_argument(
        "-C",
        "--components",
        type=int,
        choices=(1, 2, 3),
        required=True,
        help="number of galaxy components: 1=E, 2=S0/S, or 3=barred galaxy",
    )
    parser.add_argument(
        "-zp",
        "--zeropoint",
        type=float,
        default=22.5,
        help="photometric zero point (default: 22.5)",
    )
    parser.add_argument(
        "-sk",
        "--sky",
        type=float,
        default=0.0,
        help="sky background to subtract (default: 0)",
    )
    parser.add_argument(
        "-p",
        "--plate",
        type=positive_float,
        default=0.262,
        help="plate scale in arcsec/pixel (default: 0.262)",
    )
    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="take the center from the DS9 region",
    )
    parser.add_argument(
        "-no",
        "--noprint",
        action="store_true",
        help="do not print the Sérsic components from getSersic",
    )
    parser.add_argument("-m", "--mask", help="mask FITS file")
    parser.add_argument(
        "-o",
        "--out",
        default="sersic.init",
        help=(
            "base output name; an index is inserted before the extension "
            "(default: sersic.init -> sersic_001.init, ...)"
        ),
    )
    parser.add_argument(
        "-b",
        "--bards9",
        help=(
            "existing DS9 bar ellipse for the three-component model; when "
            "omitted, two half-size bar ellipses are derived from RegFile"
        ),
    )
    parser.add_argument(
        "-g",
        "--galfit_out",
        help="GALFIT output-image name forwarded to every getSersic call",
    )
    parser.add_argument(
        "-f",
        "--freeser",
        action="store_true",
        help="keep the Sérsic index of the bar component free",
    )
    parser.add_argument(
        "-cb",
        "--consbulge",
        action="store_true",
        help="add the bulge/bar axis-ratio constraints",
    )
    parser.add_argument(
        "--bar-scale",
        type=positive_float,
        default=0.5,
        help="scale applied to both RegFile ellipse radii (default: 0.5)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        help="random seed for reproducible parameter files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the getSersic arguments without executing the function",
    )
    return parser


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def compact_float(value: float) -> str:
    """Return a precise, compact representation accepted by argparse."""
    return f"{value:.10g}"


def indexed_output_name(base_name: str, index: int, count: int) -> Path:
    """Insert a one-based index before the final filename extension."""
    path = Path(base_name)
    width = max(3, len(str(count)))
    suffix = "".join(path.suffixes)
    stem = path.name[: -len(suffix)] if suffix else path.name
    return path.with_name(f"{stem}_{index:0{width}d}{suffix}")


def read_image_ellipse(region_file: Path) -> tuple[str, str, str, float, float, float]:
    """Read the first ellipse in a physical/image DS9 region file."""
    text = region_file.read_text(encoding="utf-8")
    active_coordinates: str | None = None

    for line in text.splitlines():
        coordinate_matches = list(COORDINATE_PATTERN.finditer(line))
        ellipse_match = ELLIPSE_PATTERN.search(line)

        if ellipse_match is None:
            if coordinate_matches:
                active_coordinates = coordinate_matches[-1].group(1).lower()
            continue

        coordinates = active_coordinates
        for coordinate_match in coordinate_matches:
            if coordinate_match.start() < ellipse_match.start():
                coordinates = coordinate_match.group(1).lower()

        if coordinates not in {"physical", "image"}:
            found = coordinates if coordinates is not None else "unspecified"
            raise ValueError(
                f"the ellipse coordinate system must be physical or image; found {found}"
            )

        rx = float(ellipse_match.group("rx"))
        ry = float(ellipse_match.group("ry"))
        angle = float(ellipse_match.group("angle"))
        if rx <= 0 or ry <= 0:
            raise ValueError("the ellipse radii must be greater than zero")

        return (
            coordinates,
            ellipse_match.group("x"),
            ellipse_match.group("y"),
            rx,
            ry,
            angle,
        )

    raise ValueError(f"no DS9 ellipse was found in {region_file}")


def make_bar_region_texts(region_file: Path, scale: float) -> tuple[str, str]:
    """Create half-size bar ellipses at the input PA and PA + 90 degrees."""
    coordinates, x, y, rx, ry, angle = read_image_ellipse(region_file)
    scaled_rx = compact_float(rx * scale)
    scaled_ry = compact_float(ry * scale)
    original_angle = compact_float(angle % 180.0)
    rotated_angle = compact_float((angle + 90.0) % 180.0)
    header = "# Region file format: DS9 version 4.1\n"

    original = (
        f"{header}{coordinates}\n"
        f"ellipse({x},{y},{scaled_rx},{scaled_ry},{original_angle})\n"
    )
    rotated = (
        f"{header}{coordinates}\n"
        f"ellipse({x},{y},{scaled_rx},{scaled_ry},{rotated_angle})\n"
    )
    return original, rotated


def draw_random_model_parameters(
    components: int,
    rng: random.Random,
) -> tuple[float, float | None, float | None]:
    """Draw the model-specific parameters passed to getSersic."""
    if components == 1:
        nser = rng.uniform(1.0, 6.0)
        bulgetot = None
        bulgebarat = None
    elif components == 2:
        nser = rng.uniform(1.5, 4.0)
        bulgetot = rng.uniform(0.2, 0.7)
        bulgebarat = None
    else:
        nser = rng.uniform(1.0, 3.0)
        bulgetot = rng.uniform(0.2, 0.7)
        bulgebarat = rng.uniform(0.5, 1.0)

    return nser, bulgetot, bulgebarat


def get_getsersic_arguments(
    args: argparse.Namespace,
    output: Path,
    nser: float,
    bulgetot: float | None,
    bulgebarat: float | None,
    bar_region: Path | None,
) -> list[object]:
    """Return arguments in the order required by the getSersic function."""
    return [
        args.Image,
        args.RegFile,
        args.center,
        args.mask,
        args.zeropoint,
        args.sky,
        args.noprint,
        bulgetot,
        str(bar_region) if bar_region is not None else None,
        args.plate,
        str(output),
        args.galfit_out,
        args.freeser,
        args.consbulge,
        nser,
        bulgebarat if bulgebarat is not None else 1.0,
    ]


def validate_arguments(
    parser: argparse.ArgumentParser, args: argparse.Namespace
) -> None:
    required_files = [("Image", Path(args.Image)), ("RegFile", Path(args.RegFile))]
    if args.mask:
        required_files.append(("mask", Path(args.mask)))
    if args.bards9:
        required_files.append(("bards9", Path(args.bards9)))

    for label, path in required_files:
        if not path.is_file():
            parser.error(f"{label} file does not exist: {path}")

    if args.components != 3 and args.bards9:
        parser.error("--bards9 can only be used with --components 3")
    if args.components != 3 and args.freeser:
        parser.error("--freeser can only be used with --components 3")
    if args.components != 3 and args.bar_scale != 0.5:
        parser.error("--bar-scale can only be changed with --components 3")

    output_parent = Path(args.out).expanduser().parent
    if not output_parent.is_dir():
        parser.error(f"output directory does not exist: {output_parent}")


def print_model_summary(
    index: int,
    output: Path,
    nser: float,
    bulgetot: float | None,
    bulgebarat: float | None,
    orientation: str | None,
) -> None:
    fields = [f"n={nser:.6f}"]
    if bulgetot is not None:
        fields.append(f"B/T={bulgetot:.6f}")
    if bulgebarat is not None:
        fields.append(f"bulge/bar={bulgebarat:.6f}")
    if orientation is not None:
        fields.append(f"bar={orientation}")
    print(f"[{index:03d}] {output}: " + ", ".join(fields))


def load_getsersic_runner() -> Callable[..., object]:
    """Load the core GALFITools getSersic function."""
    try:
        from galfitools.galin.getSersic import getSersic
    except ImportError as error:
        raise RuntimeError(
            "GALFITools is not installed or getSersic cannot be imported"
        ) from error
    return getSersic


def run(
    args: argparse.Namespace,
    getsersic_runner: Callable[..., object] | None = None,
) -> int:
    rng = random.Random(args.seed)
    explicit_bar = Path(args.bards9) if args.bards9 else None
    if not args.dry_run and getsersic_runner is None:
        getsersic_runner = load_getsersic_runner()

    with tempfile.TemporaryDirectory(prefix="random_getsersic_") as temp_name:
        bar_regions: tuple[Path, Path] | None = None
        if args.components == 3 and explicit_bar is None:
            original_text, rotated_text = make_bar_region_texts(
                Path(args.RegFile), args.bar_scale
            )
            temp_dir = Path(temp_name)
            original_path = temp_dir / "bar_original.reg"
            rotated_path = temp_dir / "bar_rotated_90.reg"
            original_path.write_text(original_text, encoding="utf-8")
            rotated_path.write_text(rotated_text, encoding="utf-8")
            bar_regions = (original_path, rotated_path)

        for index in range(1, args.number + 1):
            output = indexed_output_name(args.out, index, args.number)
            orientation: str | None = None

            if explicit_bar is not None:
                bar_region = explicit_bar
                orientation = "provided"
            elif bar_regions is not None:
                selected = rng.randrange(2)
                bar_region = bar_regions[selected]
                orientation = "original" if selected == 0 else "rotated90"
            else:
                bar_region = None

            nser, bulgetot, bulgebarat = draw_random_model_parameters(
                args.components, rng
            )
            function_arguments = get_getsersic_arguments(
                args,
                output,
                nser,
                bulgetot,
                bulgebarat,
                bar_region,
            )
            print_model_summary(
                index,
                output,
                nser,
                bulgetot,
                bulgebarat,
                orientation,
            )

            if args.dry_run:
                print("  getSersic(*" + repr(function_arguments) + ")")
                continue

            if getsersic_runner is None:
                raise RuntimeError("getSersic function is unavailable")
            getsersic_runner(*function_arguments)

    return 0


def mainRandomSersic(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_arguments(parser, args)

    try:
        return run(args)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(mainRandomSersic())
