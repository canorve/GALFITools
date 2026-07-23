#!/usr/bin/env python3

"""Create GALFIT mock images using sky blocks and an optional sigma image.

The GALFIT output file is assumed to contain, by default:

    HDU 1: original galaxy image
    HDU 2: GALFIT model image

Sky-only pixels are selected with DS9 ``physical`` or ``image`` box regions.
The script constructs a full-size correlated sky realization by resampling
square blocks from those regions and adds it to the GALFIT model.

When ``--sigma-image`` is supplied, the sigma map is interpreted as the total
per-pixel uncertainty. The script estimates the sky RMS from the selected sky
boxes and adds only the remaining variance:

    sigma_extra^2 = max(sigma_total^2 - sigma_sky^2, 0)

Thus, the mock image is approximately:

    mock = model + correlated_sky + Gaussian(0, sigma_extra)

Without ``--sigma-image``, the behavior is unchanged:

    mock = model + correlated_sky
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip


PROGRAM_VERSION = "3.2-physical-sigma-safe"
DEFAULT_IMAGE_EXTENSION = 1
DEFAULT_MODEL_EXTENSION = 2


def clean_extension_header(header: fits.Header) -> fits.Header:
    """Return a primary-HDU-compatible copy of an extension header."""
    output_header = header.copy()

    for keyword in (
        "XTENSION",
        "PCOUNT",
        "GCOUNT",
        "EXTNAME",
        "EXTVER",
        "EXTLEVEL",
    ):
        output_header.remove(keyword, ignore_missing=True, remove_all=True)

    return output_header


def physical_to_image_coordinates(
    points: np.ndarray,
    header: fits.Header,
) -> np.ndarray:
    """Convert DS9 physical coordinates to FITS image coordinates.

    DS9 physical coordinates use the IRAF LTM/LTV linear transformation:

        image = LTM @ physical + LTV

    The returned coordinates remain in the one-based FITS/DS9 convention.
    """
    transform = np.array(
        [
            [
                float(header.get("LTM1_1", 1.0)),
                float(header.get("LTM1_2", 0.0)),
            ],
            [
                float(header.get("LTM2_1", 0.0)),
                float(header.get("LTM2_2", 1.0)),
            ],
        ],
        dtype=float,
    )
    translation = np.array(
        [
            float(header.get("LTV1", 0.0)),
            float(header.get("LTV2", 0.0)),
        ],
        dtype=float,
    )

    return points @ transform.T + translation


def box_vertices(
    x_center: float,
    y_center: float,
    width: float,
    height: float,
    angle_degrees: float,
) -> np.ndarray:
    """Return the four vertices of a DS9 box."""
    if width <= 0 or height <= 0:
        raise ValueError("DS9 box width and height must be positive.")

    half_width = width / 2.0
    half_height = height / 2.0
    corners = np.array(
        [
            [-half_width, -half_height],
            [half_width, -half_height],
            [half_width, half_height],
            [-half_width, half_height],
        ],
        dtype=float,
    )

    angle = np.deg2rad(angle_degrees)
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle)],
            [np.sin(angle), np.cos(angle)],
        ],
        dtype=float,
    )

    return corners @ rotation.T + np.array([x_center, y_center])


def parse_ds9_box_regions(
    region_file: Path,
) -> list[tuple[bool, str, np.ndarray]]:
    """Read DS9 box regions in physical or image coordinates.

    A minus sign before ``box`` marks an excluded region. Examples::

        physical
        box(108,383.5,38,41,0)
        -box(108,383.5,5,5,0)
    """
    coordinate_systems = {
        "physical",
        "image",
        "fk4",
        "fk5",
        "icrs",
        "galactic",
        "ecliptic",
        "linear",
    }
    box_pattern = re.compile(
        r"^\s*([+-]?)\s*box\s*\(\s*([^)]*)\)\s*$",
        flags=re.IGNORECASE,
    )

    current_system = "physical"
    parsed_regions: list[tuple[bool, str, np.ndarray]] = []

    for line_number, raw_line in enumerate(
        region_file.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        # DS9 properties follow '#'. They are not needed for the mask.
        line = line.split("#", maxsplit=1)[0].strip()

        if not line or line.lower().startswith("global"):
            continue

        for statement in line.split(";"):
            statement = statement.strip()

            if not statement:
                continue

            lower_statement = statement.lower()

            if lower_statement in coordinate_systems:
                current_system = lower_statement
                continue

            match = box_pattern.match(statement)

            if match is None:
                raise ValueError(
                    f"Unsupported DS9 region at line {line_number}: "
                    f"{statement!r}. This script accepts box regions only."
                )

            if current_system not in {"physical", "image"}:
                raise ValueError(
                    f"Unsupported coordinate system {current_system!r} at "
                    f"line {line_number}. Use DS9 physical or image boxes."
                )

            values = [value.strip() for value in match.group(2).split(",")]

            if len(values) != 5:
                raise ValueError(
                    f"Invalid DS9 box at line {line_number}: expected five "
                    "values x, y, width, height, angle."
                )

            try:
                x_center, y_center, width, height, angle = map(float, values)
            except ValueError as error:
                raise ValueError(
                    f"Invalid numeric value in DS9 box at line {line_number}."
                ) from error

            include = match.group(1) != "-"
            vertices = box_vertices(
                x_center=x_center,
                y_center=y_center,
                width=width,
                height=height,
                angle_degrees=angle,
            )
            parsed_regions.append((include, current_system, vertices))

    if not parsed_regions:
        raise ValueError(f"No DS9 box regions were found in {region_file}.")

    return parsed_regions


def polygon_to_mask(
    vertices: np.ndarray,
    shape: tuple[int, int],
) -> np.ndarray:
    """Rasterize a convex polygon using NumPy pixel-center coordinates."""
    ny, nx = shape
    x_min = max(0, int(np.floor(np.min(vertices[:, 0]))))
    x_max = min(nx - 1, int(np.ceil(np.max(vertices[:, 0]))))
    y_min = max(0, int(np.floor(np.min(vertices[:, 1]))))
    y_max = min(ny - 1, int(np.ceil(np.max(vertices[:, 1]))))

    output = np.zeros(shape, dtype=bool)

    if x_min > x_max or y_min > y_max:
        return output

    x_grid, y_grid = np.meshgrid(
        np.arange(x_min, x_max + 1, dtype=float),
        np.arange(y_min, y_max + 1, dtype=float),
    )

    cross_products = []

    for index in range(len(vertices)):
        x0, y0 = vertices[index]
        x1, y1 = vertices[(index + 1) % len(vertices)]
        cross_products.append((x1 - x0) * (y_grid - y0) - (y1 - y0) * (x_grid - x0))

    cross_products = np.asarray(cross_products)
    tolerance = 1.0e-10
    inside = np.all(cross_products >= -tolerance, axis=0) | np.all(
        cross_products <= tolerance,
        axis=0,
    )

    output[y_min : y_max + 1, x_min : x_max + 1] = inside
    return output


def read_region_mask(
    region_file: Path,
    shape: tuple[int, int],
    header: fits.Header,
) -> np.ndarray:
    """Create a Boolean mask from DS9 physical or image box regions."""
    included_mask = np.zeros(shape, dtype=bool)
    excluded_mask = np.zeros(shape, dtype=bool)
    used_regions = 0

    for include, coordinate_system, vertices in parse_ds9_box_regions(region_file):
        if coordinate_system == "physical":
            image_vertices = physical_to_image_coordinates(vertices, header)
        else:
            image_vertices = vertices

        # DS9/FITS pixel centers are one-based; NumPy pixel centers are zero-based.
        pixel_vertices = image_vertices - 1.0
        current_mask = polygon_to_mask(pixel_vertices, shape)

        if not np.any(current_mask):
            continue

        if include:
            included_mask |= current_mask
        else:
            excluded_mask |= current_mask

        used_regions += 1

    combined_mask = included_mask & ~excluded_mask

    if used_regions == 0 or not np.any(combined_mask):
        raise ValueError(
            "None of the included DS9 box regions overlap the input image."
        )

    return combined_mask


def estimate_sky_level(
    image: np.ndarray,
    sky_mask: np.ndarray,
    sigma: float = 3.0,
    maxiters: int = 5,
) -> tuple[float, float, int]:
    """Estimate the sky level and RMS from sigma-clipped sky pixels."""
    sky_pixels = image[sky_mask & np.isfinite(image)]

    if sky_pixels.size == 0:
        raise ValueError("No finite sky pixels are available.")

    clipped = sigma_clip(
        sky_pixels,
        sigma=sigma,
        maxiters=maxiters,
        masked=True,
    )
    good_pixels = np.asarray(clipped.compressed(), dtype=float)

    if good_pixels.size < 2:
        raise ValueError("Too few sky pixels remain after sigma clipping.")

    sky_level = float(np.median(good_pixels))
    sky_rms = float(np.std(good_pixels, ddof=1))

    if not np.isfinite(sky_rms) or sky_rms <= 0:
        raise ValueError("The estimated sky RMS is not positive and finite.")

    return sky_level, sky_rms, int(good_pixels.size)


def load_2d_fits_image(
    filename: Path,
    extension: int | None = None,
) -> tuple[np.ndarray, fits.Header, int]:
    """Load a 2D FITS image from a specified HDU or the first 2D HDU."""
    filename = filename.expanduser().resolve()

    if not filename.is_file():
        raise FileNotFoundError(f"FITS image not found: {filename}")

    with fits.open(filename, memmap=False) as hdulist:
        if extension is not None:
            if extension < 0 or extension >= len(hdulist):
                raise ValueError(
                    f"FITS file {filename} does not contain HDU {extension}."
                )

            hdu = hdulist[extension]

            if hdu.data is None:
                raise ValueError(
                    f"HDU {extension} in {filename} contains no image data."
                )

            data = np.asarray(hdu.data, dtype=float)

            if data.ndim != 2:
                raise ValueError(f"HDU {extension} in {filename} is not a 2D image.")

            return data, hdu.header.copy(), extension

        for hdu_index, hdu in enumerate(hdulist):
            if hdu.data is None:
                continue

            data = np.asarray(hdu.data, dtype=float)

            if data.ndim == 2:
                return data, hdu.header.copy(), hdu_index

    raise ValueError(f"No 2D image was found in {filename}.")


def infer_noise_map_kind(
    filename: Path,
    header: fits.Header,
) -> str:
    """Infer whether a noise map contains sigma or inverse variance."""
    tokens = " ".join(
        [
            filename.name.lower(),
            str(header.get("EXTNAME", "")).lower(),
            str(header.get("BUNIT", "")).lower(),
            str(header.get("IMTYPE", "")).lower(),
        ]
    )

    invvar_tokens = ("invvar", "ivar", "inverse variance", "weight", "wht")
    sigma_tokens = ("sigma", "std", "stdev", "rms", "error", "err")

    if any(token in tokens for token in invvar_tokens):
        return "invvar"

    if any(token in tokens for token in sigma_tokens):
        return "sigma"

    print(
        "Warning: could not infer the noise-map type from its filename or "
        "header. Assuming sigma. Use --sigma-kind explicitly if necessary."
    )
    return "sigma"


def convert_noise_map_to_sigma(
    data: np.ndarray,
    kind: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a sigma or inverse-variance image into a sigma map.

    Returns
    -------
    sigma_map
        Sigma values. Invalid pixels are set to zero.
    valid_mask
        True where the original noise-map value was valid.
    """
    data = np.asarray(data, dtype=float)

    if kind == "sigma":
        valid_mask = np.isfinite(data) & (data > 0)
        sigma_map = np.zeros_like(data, dtype=float)
        sigma_map[valid_mask] = data[valid_mask]
        return sigma_map, valid_mask

    if kind == "invvar":
        valid_mask = np.isfinite(data) & (data > 0)
        sigma_map = np.zeros_like(data, dtype=float)
        sigma_map[valid_mask] = 1.0 / np.sqrt(data[valid_mask])
        return sigma_map, valid_mask

    raise ValueError(f"Unsupported noise-map kind: {kind!r}.")


def compute_extra_sigma_map(
    total_sigma: np.ndarray,
    valid_sigma: np.ndarray,
    sky_rms: float,
) -> tuple[np.ndarray, int, int]:
    """Compute the non-sky sigma after subtracting sky variance.

    The calculation is

        sigma_extra = sqrt(max(sigma_total**2 - sky_rms**2, 0)).

    Invalid sigma-map pixels receive zero extra noise.
    """
    extra_sigma = np.zeros_like(total_sigma, dtype=float)
    total_variance = total_sigma[valid_sigma] ** 2
    extra_variance = total_variance - sky_rms**2

    clipped_count = int(np.count_nonzero(extra_variance <= 0))
    extra_sigma[valid_sigma] = np.sqrt(np.clip(extra_variance, 0.0, None))
    invalid_count = int(np.count_nonzero(~valid_sigma))

    return extra_sigma, clipped_count, invalid_count


def find_valid_block_origins(
    valid_mask: np.ndarray,
    block_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Find top-left positions of square blocks fully inside valid pixels."""
    ny, nx = valid_mask.shape

    if block_size > ny or block_size > nx:
        return np.array([], dtype=int), np.array([], dtype=int)

    integral = np.pad(
        valid_mask.astype(np.int32),
        ((1, 0), (1, 0)),
        mode="constant",
    )
    integral = integral.cumsum(axis=0).cumsum(axis=1)

    block_sums = (
        integral[block_size:, block_size:]
        - integral[:-block_size, block_size:]
        - integral[block_size:, :-block_size]
        + integral[:-block_size, :-block_size]
    )

    return np.where(block_sums == block_size * block_size)


def select_block_size(
    valid_mask: np.ndarray,
    requested_size: int,
) -> tuple[int, np.ndarray, np.ndarray]:
    """Select the largest usable block size not exceeding the request."""
    if requested_size < 2:
        raise ValueError("The block size must be at least 2 pixels.")

    maximum_size = min(requested_size, *valid_mask.shape)

    for block_size in range(maximum_size, 1, -1):
        y_origins, x_origins = find_valid_block_origins(
            valid_mask,
            block_size,
        )

        if y_origins.size > 0:
            return block_size, y_origins, x_origins

    raise ValueError("The sky regions do not contain a contiguous 2 x 2 pixel block.")


def transform_square_block(
    block: np.ndarray,
    rotation: int,
    reflect: bool,
) -> np.ndarray:
    """Apply a specified rotation and optional reflection to a square block."""
    transformed = np.rot90(block, k=rotation)

    if reflect:
        transformed = np.fliplr(transformed)

    return transformed


def create_sky_realization(
    noise_source: np.ndarray,
    output_shape: tuple[int, int],
    block_size: int,
    y_origins: np.ndarray,
    x_origins: np.ndarray,
    rng: np.random.Generator,
    transform_blocks: bool = False,
) -> np.ndarray:
    """Fill an image using randomly sampled correlated-noise blocks."""
    ny, nx = output_shape
    sky_image = np.empty(output_shape, dtype=float)

    for y_out in range(0, ny, block_size):
        height = min(block_size, ny - y_out)

        for x_out in range(0, nx, block_size):
            width = min(block_size, nx - x_out)
            choice = int(rng.integers(0, y_origins.size))
            y_in = int(y_origins[choice])
            x_in = int(x_origins[choice])

            block = noise_source[
                y_in : y_in + block_size,
                x_in : x_in + block_size,
            ].copy()

            if transform_blocks:
                rotation = int(rng.integers(0, 4))
                reflect = bool(rng.integers(0, 2))
                block = transform_square_block(block, rotation, reflect)

            sky_image[
                y_out : y_out + height,
                x_out : x_out + width,
            ] = block[:height, :width]

    return sky_image


def create_mock_images(
    galfit_file: Path,
    region_file: Path,
    output_directory: Path,
    number: int = 100,
    block_size: int = 16,
    image_extension: int = DEFAULT_IMAGE_EXTENSION,
    model_extension: int = DEFAULT_MODEL_EXTENSION,
    prefix: str | None = None,
    seed: int | None = None,
    keep_sky_level: bool = False,
    transform_blocks: bool = False,
    save_sky: bool = True,
    sigma_image: Path | None = None,
    sigma_extension: int | None = None,
    sigma_kind: str = "auto",
    sigma_max_factor: float = 100.0,
) -> list[Path]:
    """Create mock images from a GALFIT model and sampled sky noise."""
    galfit_file = galfit_file.expanduser().resolve()
    region_file = region_file.expanduser().resolve()
    output_directory = output_directory.expanduser().resolve()

    if not galfit_file.is_file():
        raise FileNotFoundError(f"GALFIT file not found: {galfit_file}")

    if not region_file.is_file():
        raise FileNotFoundError(f"Region file not found: {region_file}")

    if number < 1:
        raise ValueError("The number of mock images must be at least 1.")

    output_directory.mkdir(parents=True, exist_ok=True)

    if prefix is None:
        prefix = galfit_file.stem

    with fits.open(galfit_file, memmap=False) as hdulist:
        required_extension = max(image_extension, model_extension)

        if len(hdulist) <= required_extension:
            raise ValueError(
                f"The GALFIT file does not contain HDU {required_extension}."
            )

        image_hdu = hdulist[image_extension]
        model_hdu = hdulist[model_extension]

        if image_hdu.data is None:
            raise ValueError(f"HDU {image_extension} does not contain image data.")

        if model_hdu.data is None:
            raise ValueError(f"HDU {model_extension} does not contain model data.")

        image = np.asarray(image_hdu.data, dtype=float)
        model = np.asarray(model_hdu.data, dtype=float)
        image_header = image_hdu.header.copy()
        output_header = clean_extension_header(model_hdu.header)

    if image.ndim != 2 or model.ndim != 2:
        raise ValueError("The original and model images must be two-dimensional.")

    if image.shape != model.shape:
        raise ValueError(
            "The original and model images have different shapes: "
            f"{image.shape} and {model.shape}."
        )

    sky_mask = read_region_mask(
        region_file=region_file,
        shape=image.shape,
        header=image_header,
    )
    sky_level, sky_rms, clipped_sky_pixels = estimate_sky_level(
        image,
        sky_mask,
    )

    # The donor field contains sky fluctuations around zero. This prevents
    # adding the mean sky twice when the GALFIT model includes a sky component.
    noise_source = image - sky_level
    valid_block_mask = sky_mask & np.isfinite(noise_source)

    total_sigma: np.ndarray | None = None
    extra_sigma: np.ndarray | None = None
    sigma_file: Path | None = None
    used_sigma_extension: int | None = None
    clipped_variance_pixels = 0
    invalid_sigma_pixels = 0
    extreme_sigma_pixels = 0
    effective_sigma_kind = sigma_kind
    median_sigma_sky = 0.0
    sigma_limit = 0.0

    if sigma_image is not None:
        sigma_file = sigma_image.expanduser().resolve()
        sigma_data, sigma_header, used_sigma_extension = load_2d_fits_image(
            sigma_file,
            sigma_extension,
        )

        if sigma_data.shape != image.shape:
            raise ValueError(
                "The sigma/noise image shape does not match the GALFIT "
                f"images: {sigma_data.shape} and {image.shape}."
            )

        effective_sigma_kind = sigma_kind
        if effective_sigma_kind == "auto":
            effective_sigma_kind = infer_noise_map_kind(
                sigma_file,
                sigma_header,
            )

        total_sigma, valid_sigma = convert_noise_map_to_sigma(
            sigma_data,
            effective_sigma_kind,
        )

        valid_sky_sigma = total_sigma[sky_mask & valid_sigma]
        if valid_sky_sigma.size == 0:
            raise ValueError(
                "The sigma image has no valid pixels inside the DS9 sky boxes."
            )

        median_sigma_sky = float(np.median(valid_sky_sigma))
        if not np.isfinite(median_sigma_sky) or median_sigma_sky <= 0:
            raise ValueError("The median sigma inside the DS9 sky boxes is invalid.")

        # Tiny positive inverse-variance values can produce enormous sigma
        # values in poorly covered pixels. They dominate DS9 scaling and can
        # hide the model. Treat such pixels as invalid for extra-noise
        # generation. This threshold is intentionally generous.
        if sigma_max_factor <= 1:
            raise ValueError("--sigma-max-factor must be greater than 1.")

        sigma_limit = sigma_max_factor * median_sigma_sky
        extreme_sigma = valid_sigma & (total_sigma > sigma_limit)
        extreme_sigma_pixels = int(np.count_nonzero(extreme_sigma))
        valid_sigma[extreme_sigma] = False
        total_sigma[extreme_sigma] = 0.0

        # The block-resampled sky already contributes the measured sky RMS.
        # Add only the remaining variance from the total sigma map.
        (
            extra_sigma,
            clipped_variance_pixels,
            invalid_sigma_pixels,
        ) = compute_extra_sigma_map(
            total_sigma=total_sigma,
            valid_sigma=valid_sigma,
            sky_rms=sky_rms,
        )

        # Extra source-related variance should not be generated in invalid
        # model pixels.
        extra_sigma[~np.isfinite(model)] = 0.0

        positive_sigma = total_sigma[total_sigma > 0]
        if positive_sigma.size:
            median_total_sigma = float(np.median(positive_sigma))
            if median_total_sigma > 20.0 * max(sky_rms, median_sigma_sky):
                raise ValueError(
                    "The noise map is much larger than the measured sky "
                    "noise. It may be an inverse-variance map being read as "
                    "sigma. Retry with --sigma-kind invvar."
                )

    used_block_size, y_origins, x_origins = select_block_size(
        valid_block_mask,
        block_size,
    )

    if used_block_size != block_size:
        print(
            f"Warning: requested block size {block_size} is not available. "
            f"Using {used_block_size} pixels instead."
        )

    print(f"Sky-region pixels:       {int(np.count_nonzero(sky_mask))}")
    print(f"Clipped sky pixels:      {clipped_sky_pixels}")
    print(f"Sigma-clipped sky level: {sky_level:.8g}")
    print(f"Sigma-clipped sky RMS:   {sky_rms:.8g}")
    print(f"Valid donor blocks:      {y_origins.size}")
    print(f"Block size used:         {used_block_size}")

    if sigma_file is None:
        print("Sigma image:             not supplied")
        print("Mock formula:            model + correlated sky")
    else:
        valid_total = total_sigma[total_sigma > 0] if total_sigma is not None else []
        valid_extra = extra_sigma[extra_sigma > 0] if extra_sigma is not None else []
        median_total = float(np.median(valid_total)) if len(valid_total) else 0.0
        median_extra = float(np.median(valid_extra)) if len(valid_extra) else 0.0
        valid_count = int(np.count_nonzero(total_sigma > 0))
        clipped_fraction = (
            100.0 * clipped_variance_pixels / valid_count if valid_count else 0.0
        )

        print(f"Sigma image:             {sigma_file}")
        print(f"Sigma HDU:               {used_sigma_extension}")
        print(f"Sigma input kind:        {effective_sigma_kind}")
        print(f"Median sigma in sky:     {median_sigma_sky:.8g}")
        print(f"Maximum accepted sigma:  {sigma_limit:.8g}")
        print(f"Extreme sigma pixels:    {extreme_sigma_pixels}")
        print(f"Median total sigma:      {median_total:.8g}")
        print(f"Median extra sigma:      {median_extra:.8g}")
        print(f"Invalid sigma pixels:    {invalid_sigma_pixels}")
        print(
            "Pixels with no extra variance: "
            f"{clipped_variance_pixels} ({clipped_fraction:.2f}%)"
        )
        print("Mock formula:            model + correlated sky + extra noise")

    rng = np.random.default_rng(seed)
    output_files: list[Path] = []

    for index in range(1, number + 1):
        sky_noise = create_sky_realization(
            noise_source=noise_source,
            output_shape=image.shape,
            block_size=used_block_size,
            y_origins=y_origins,
            x_origins=x_origins,
            rng=rng,
            transform_blocks=transform_blocks,
        )

        if keep_sky_level:
            sky_image = sky_noise + sky_level
        else:
            sky_image = sky_noise

        if extra_sigma is None:
            extra_noise = np.zeros_like(model, dtype=float)
        else:
            extra_noise = (
                rng.normal(
                    loc=0.0,
                    scale=1.0,
                    size=model.shape,
                )
                * extra_sigma
            )

        mock_image = model + sky_image + extra_noise

        if index == 1:

            def summarize(name: str, array: np.ndarray) -> None:
                finite = np.asarray(array, dtype=float)
                finite = finite[np.isfinite(finite)]
                if finite.size == 0:
                    print(f"{name:24s}: no finite pixels")
                    return
                p01, p50, p99 = np.percentile(finite, [1, 50, 99])
                print(
                    f"{name:24s}: min={finite.min():.6g}, "
                    f"p01={p01:.6g}, median={p50:.6g}, "
                    f"p99={p99:.6g}, max={finite.max():.6g}"
                )

            print("First-realization diagnostics:")
            summarize("GALFIT model", model)
            summarize("Sky realization", sky_image)
            summarize("Extra Gaussian noise", extra_noise)
            summarize("Final mock", mock_image)

        header = output_header.copy()
        header["MOCKNUM"] = (index, "Mock realization number")
        header["BLKSIZE"] = (used_block_size, "Sky bootstrap block size")
        header["SKYLEVEL"] = (sky_level, "Median sky in selected regions")
        header["SKYRMS"] = (sky_rms, "Sigma-clipped sky RMS")
        header["KEEPSKY"] = (keep_sky_level, "Absolute sky level retained")
        header["IMGEXT"] = (image_extension, "Original-image HDU")
        header["MODELEXT"] = (model_extension, "GALFIT model HDU")
        header["RNGSEED"] = (-1 if seed is None else seed, "Base RNG seed")
        header["SIGUSED"] = (sigma_file is not None, "Sigma/noise map used")

        if sigma_file is not None:
            header["SIGFILE"] = (sigma_file.name, "Sigma/noise map filename")
            header["SIGEXT"] = (
                int(used_sigma_extension),
                "Sigma/noise map HDU",
            )
            header["SIGKIND"] = (
                effective_sigma_kind,
                "Input map: sigma or invvar",
            )
            header["SIGMAXF"] = (
                sigma_max_factor,
                "Maximum sigma / median sky sigma",
            )
            header["SIGBAD"] = (
                extreme_sigma_pixels,
                "Extreme sigma pixels excluded",
            )
            header["VARCLIP"] = (
                clipped_variance_pixels,
                "Pixels where extra variance clipped to zero",
            )
            header.add_history(
                "Extra variance = max(total sigma^2 - measured sky RMS^2, 0)."
            )

        header.add_history(
            "Mock = GALFIT model + block-resampled correlated sky noise."
        )

        if sigma_file is not None:
            header.add_history(
                "Independent Gaussian extra noise was added from sigma map."
            )

        header.add_history(f"Sky regions: {region_file.name}")
        header.add_history(f"Source GALFIT file: {galfit_file.name}")

        if save_sky:
            sky_file = output_directory / f"{prefix}_sky_{index:04d}.fits"
            fits.PrimaryHDU(data=sky_image, header=header.copy(),).writeto(
                sky_file,
                overwrite=True,
                output_verify="fix",
            )
            output_files.append(sky_file)

        mock_file = output_directory / f"{prefix}_mock_{index:04d}.fits"
        fits.PrimaryHDU(data=mock_image, header=header,).writeto(
            mock_file,
            overwrite=True,
            output_verify="fix",
        )
        output_files.append(mock_file)

        print(f"Created realization {index:04d}: {mock_file}")

    return output_files


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Create GALFIT mock galaxy images by adding block-resampled "
            "correlated sky noise from DS9 physical/image boxes to the "
            "model. An optional sigma or inverse-variance image can add "
            "only the variance not already represented by the sky blocks."
        )
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {PROGRAM_VERSION}",
    )
    parser.add_argument(
        "galfit_file",
        type=Path,
        help="GALFIT output cube, such as img-out.fits.",
    )
    parser.add_argument(
        "region_file",
        type=Path,
        help=(
            "DS9 file containing source-free sky boxes in physical or "
            "image coordinates."
        ),
    )
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        default=100,
        help="Number of mock images. Default: 100.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("sky_mocks"),
        help="Output directory. Default: sky_mocks.",
    )
    parser.add_argument(
        "-b",
        "--block-size",
        type=int,
        default=16,
        help=(
            "Square sky-bootstrap block size in pixels. It should exceed "
            "the noise-correlation length. Default: 16."
        ),
    )
    parser.add_argument(
        "--image-ext",
        type=int,
        default=DEFAULT_IMAGE_EXTENSION,
        help="HDU containing the original galaxy image. Default: 1.",
    )
    parser.add_argument(
        "--model-ext",
        type=int,
        default=DEFAULT_MODEL_EXTENSION,
        help="HDU containing the GALFIT model image. Default: 2.",
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
        help="Random seed for reproducible realizations.",
    )
    parser.add_argument(
        "--keep-sky-level",
        action="store_true",
        help=(
            "Retain the absolute sky level in the sampled sky image. Use "
            "only when the GALFIT model does not include the sky level."
        ),
    )
    parser.add_argument(
        "--transform-blocks",
        action="store_true",
        help=(
            "Randomly rotate and reflect sampled sky blocks. Their original "
            "orientation is preserved by default."
        ),
    )
    parser.add_argument(
        "--mock-only",
        action="store_true",
        help="Write only mock images, not the intermediate sky images.",
    )
    parser.add_argument(
        "--sigma-image",
        type=Path,
        default=None,
        help=(
            "Optional GALFIT sigma image or inverse-variance image. When "
            "provided, only variance beyond the measured sky RMS is added."
        ),
    )
    parser.add_argument(
        "--sigma-ext",
        type=int,
        default=None,
        help=(
            "HDU containing the sigma/noise map. By default, the first 2D "
            "image HDU is used."
        ),
    )
    parser.add_argument(
        "--sigma-kind",
        choices=("auto", "sigma", "invvar"),
        default="auto",
        help=(
            "Interpret --sigma-image as sigma or inverse variance. "
            "Default: auto, inferred from filename/header."
        ),
    )
    parser.add_argument(
        "--sigma-max-factor",
        type=float,
        default=100.0,
        help=(
            "Ignore sigma values larger than this factor times the median "
            "sigma measured inside the sky boxes. This suppresses poorly "
            "covered pixels with tiny positive inverse variance. Default: 100."
        ),
    )

    return parser.parse_args()


def mainGalfitSkyMock() -> None:
    """Run the command-line program."""
    args = parse_arguments()
    # print(f"create_galfit_sky_mocks version {PROGRAM_VERSION}")
    print("DS9 parser: built-in physical/image box parser (no regions package)")

    try:
        create_mock_images(
            galfit_file=args.galfit_file,
            region_file=args.region_file,
            output_directory=args.output_dir,
            number=args.number,
            block_size=args.block_size,
            image_extension=args.image_ext,
            model_extension=args.model_ext,
            prefix=args.prefix,
            seed=args.seed,
            keep_sky_level=args.keep_sky_level,
            transform_blocks=args.transform_blocks,
            save_sky=not args.mock_only,
            sigma_image=args.sigma_image,
            sigma_extension=args.sigma_ext,
            sigma_kind=args.sigma_kind,
            sigma_max_factor=args.sigma_max_factor,
        )
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error


if __name__ == "__main__":
    mainGalfitSkyMock()
