#!/usr/bin/env python3

"""Separate a GALFIT output FITS file into input, model, and residual images."""

from __future__ import annotations

import argparse
from pathlib import Path

from astropy.io import fits


GALFIT_EXTENSIONS = {
    "galaxy": 1,
    "model": 2,
    "residual": 3,
}


def extension_to_primary_hdu(
    extension: fits.ImageHDU,
) -> fits.PrimaryHDU:
    """
    Convert a FITS image extension into a standalone primary HDU.

    The image data and corresponding extension header are preserved.
    Extension-only structural keywords are removed.
    """
    header = extension.header.copy()

    extension_keywords = (
        "XTENSION",
        "PCOUNT",
        "GCOUNT",
        "EXTNAME",
        "EXTVER",
        "EXTLEVEL",
    )

    for keyword in extension_keywords:
        header.remove(keyword, ignore_missing=True, remove_all=True)

    return fits.PrimaryHDU(
        data=extension.data,
        header=header,
    )


def split_galfit_output(
    input_file: Path,
    output_directory: Path | None = None,
    prefix: str | None = None,
    overwrite: bool = False,
) -> list[Path]:
    """Separate a GALFIT output file into its three image products."""
    input_file = input_file.expanduser().resolve()

    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    if output_directory is None:
        output_directory = input_file.parent
    else:
        output_directory = output_directory.expanduser().resolve()

    output_directory.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = input_file.stem

    output_files: list[Path] = []

    with fits.open(input_file, memmap=False) as hdulist:
        if len(hdulist) < 4:
            raise ValueError(
                "The GALFIT file must contain at least four HDUs: "
                "a primary HDU and image extensions 1, 2, and 3."
            )

        for image_name, extension_number in GALFIT_EXTENSIONS.items():
            extension = hdulist[extension_number]

            if extension.data is None:
                raise ValueError(
                    f"HDU {extension_number}, expected to contain the "
                    f"{image_name} image, has no data."
                )

            output_file = output_directory / f"{prefix}_{image_name}.fits"

            output_hdu = extension_to_primary_hdu(extension)
            output_hdu.writeto(
                output_file,
                overwrite=overwrite,
                output_verify="fix",
            )

            output_files.append(output_file)

            print(f"HDU {extension_number}: " f"{image_name:8s} -> {output_file}")

    return output_files


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Separate GALFIT output extensions into galaxy, "
            "model, and residual FITS files."
        )
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="GALFIT output FITS file.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: input-file directory.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        default=None,
        help="Output filename prefix. Default: input filename stem.",
    )

    parser.add_argument(
        "-f",
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )

    return parser.parse_args()


def mainSplitCube() -> None:
    """Run the command-line program."""
    args = parse_arguments()

    try:
        split_galfit_output(
            input_file=args.input_file,
            output_directory=args.output_dir,
            prefix=args.prefix,
            overwrite=args.overwrite,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error


if __name__ == "__main__":
    mainSplitCube()
