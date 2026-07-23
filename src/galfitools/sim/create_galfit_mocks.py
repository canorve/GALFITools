#!/usr/bin/env python3

"""Create mock galaxy images from a GALFIT model and residual image."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


MODEL_EXTENSION = 2
RESIDUAL_EXTENSION = 3


def clean_extension_header(header: fits.Header) -> fits.Header:
    """Convert an image-extension header into a primary-HDU-compatible header."""
    output_header = header.copy()

    for keyword in (
        "XTENSION",
        "PCOUNT",
        "GCOUNT",
        "EXTNAME",
        "EXTVER",
        "EXTLEVEL",
    ):
        output_header.remove(
            keyword,
            ignore_missing=True,
            remove_all=True,
        )

    return output_header


def transform_residual(
    residual: np.ndarray,
    rng: np.random.Generator,
    rotate: bool = True,
    reflect: bool = True,
) -> tuple[np.ndarray, int, int, int, bool]:
    """
    Randomly transform and circularly shift a residual image.

    Circular shifts preserve the residual pixel values and much of the
    spatially correlated noise structure.
    """
    transformed = np.array(residual, dtype=float, copy=True)

    rotation = 0
    reflected = False

    if rotate:
        rotation = int(rng.integers(0, 4))
        transformed = np.rot90(transformed, k=rotation)

    if reflect and bool(rng.integers(0, 2)):
        transformed = np.fliplr(transformed)
        reflected = True

    # A 90-degree rotation changes the array shape for non-square images.
    if transformed.shape != residual.shape:
        transformed = np.resize(transformed, residual.shape)

    ny, nx = transformed.shape

    shift_y = int(rng.integers(0, ny))
    shift_x = int(rng.integers(0, nx))

    transformed = np.roll(
        transformed,
        shift=(shift_y, shift_x),
        axis=(0, 1),
    )

    return transformed, shift_y, shift_x, rotation, reflected


def create_mock_images(
    galfit_file: Path,
    output_directory: Path,
    number: int,
    prefix: str | None = None,
    seed: int | None = None,
    rotate: bool = True,
    reflect: bool = True,
) -> list[Path]:
    """Create mock galaxy images from a GALFIT output cube."""
    galfit_file = galfit_file.expanduser().resolve()
    output_directory = output_directory.expanduser().resolve()

    if not galfit_file.is_file():
        raise FileNotFoundError(f"GALFIT file not found: {galfit_file}")

    if number < 1:
        raise ValueError("The number of mock images must be at least 1.")

    output_directory.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = galfit_file.stem

    rng = np.random.default_rng(seed)

    with fits.open(galfit_file, memmap=False) as hdulist:
        if len(hdulist) <= RESIDUAL_EXTENSION:
            raise ValueError(
                "The GALFIT file does not contain the expected model and "
                "residual extensions."
            )

        model_hdu = hdulist[MODEL_EXTENSION]
        residual_hdu = hdulist[RESIDUAL_EXTENSION]

        if model_hdu.data is None:
            raise ValueError(f"HDU {MODEL_EXTENSION} does not contain model data.")

        if residual_hdu.data is None:
            raise ValueError(
                f"HDU {RESIDUAL_EXTENSION} does not contain residual data."
            )

        model = np.asarray(model_hdu.data, dtype=float)
        residual = np.asarray(residual_hdu.data, dtype=float)

        if model.ndim != 2 or residual.ndim != 2:
            raise ValueError("The model and residual must be two-dimensional.")

        if model.shape != residual.shape:
            raise ValueError(
                "The model and residual images have different dimensions: "
                f"{model.shape} and {residual.shape}."
            )

        output_header = clean_extension_header(model_hdu.header)

    output_files: list[Path] = []

    for index in range(1, number + 1):
        mock_residual, shift_y, shift_x, rotation, reflected = transform_residual(
            residual=residual,
            rng=rng,
            rotate=rotate,
            reflect=reflect,
        )

        mock_image = model + mock_residual

        header = output_header.copy()
        header["MOCKNUM"] = (index, "Mock realization number")
        header["RESIDY"] = (shift_y, "Residual circular shift in Y")
        header["RESIDX"] = (shift_x, "Residual circular shift in X")
        header["RESROT"] = (rotation * 90, "Residual rotation in degrees")
        header["RESFLIP"] = (reflected, "Residual reflected horizontally")
        header["MOCKSEED"] = (
            -1 if seed is None else seed,
            "Random-number generator seed",
        )
        header.add_history(
            "Mock image created as GALFIT model plus transformed residual."
        )
        header.add_history(f"Source GALFIT file: {galfit_file.name}")

        output_file = output_directory / f"{prefix}_mock_{index:04d}.fits"

        fits.PrimaryHDU(data=mock_image, header=header,).writeto(
            output_file,
            overwrite=True,
            output_verify="fix",
        )

        output_files.append(output_file)

        print(
            f"{output_file}: shift=({shift_y}, {shift_x}), "
            f"rotation={rotation * 90} deg, reflected={reflected}"
        )

    return output_files


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Create mock galaxy images by adding transformed GALFIT "
            "residuals to the GALFIT model."
        )
    )

    parser.add_argument(
        "galfit_file",
        type=Path,
        help="GALFIT output cube, such as img-out.fits.",
    )

    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=100,
        help="Number of mock images to create. Default: 100.",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("mock_images"),
        help="Output directory. Default: mock_images.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        default=None,
        help="Output filename prefix. Default: GALFIT filename stem.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible mock images.",
    )

    parser.add_argument(
        "--no-rotate",
        action="store_true",
        help="Do not randomly rotate the residual image.",
    )

    parser.add_argument(
        "--no-reflect",
        action="store_true",
        help="Do not randomly reflect the residual image.",
    )

    return parser.parse_args()


def mainGalfitMock() -> None:
    """Run the command-line program."""
    args = parse_arguments()

    try:
        create_mock_images(
            galfit_file=args.galfit_file,
            output_directory=args.output_dir,
            number=args.number,
            prefix=args.prefix,
            seed=args.seed,
            rotate=not args.no_rotate,
            reflect=not args.no_reflect,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error


if __name__ == "__main__":
    mainGalfitMock()
