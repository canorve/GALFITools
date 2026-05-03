#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import label


def read_first_image_hdu(filename: str) -> tuple[np.ndarray, fits.Header]:
    """
    Read the first FITS HDU that contains image data.

    Parameters
    ----------
    filename : str
        Input FITS filename.

    Returns
    -------
    data : numpy.ndarray
        Image data array.
    header : astropy.io.fits.Header
        Header of the image HDU.
    """
    with fits.open(filename) as hdul:
        for hdu in hdul:
            if hdu.data is not None:
                return np.asarray(hdu.data), hdu.header

    raise ValueError(f"No image data found in {filename}")


def build_selected_mask(
    image: np.ndarray,
    target_value: int = 4096,
    use_bit: bool = False,
    bit: int = 12,
) -> np.ndarray:
    """
    Build a boolean mask from a maskbits image.

    Parameters
    ----------
    image : numpy.ndarray
        Input maskbits image.
    target_value : int, optional
        Exact pixel value used when use_bit is False.
    use_bit : bool, optional
        If True, select pixels where the given bit is active.
    bit : int, optional
        Bit number used when use_bit is True.

    Returns
    -------
    numpy.ndarray
        Boolean selected mask.
    """
    if use_bit:
        return (image.astype(np.int64) & (1 << bit)) != 0

    return image == target_value


def select_central_component(selected_mask: np.ndarray) -> np.ndarray:
    """
    Select the connected component closest to the image center.

    Parameters
    ----------
    selected_mask : numpy.ndarray
        Boolean mask.

    Returns
    -------
    numpy.ndarray
        Boolean mask for the selected central component.
    """
    labeled_mask, num_components = label(selected_mask)

    if num_components == 0:
        raise RuntimeError("No selected pixels were found.")

    ny, nx = selected_mask.shape
    x_image_center = (nx - 1) / 2.0
    y_image_center = (ny - 1) / 2.0

    center_label = labeled_mask[int(round(y_image_center)), int(round(x_image_center))]

    if center_label > 0:
        return labeled_mask == center_label

    best_label = None
    best_distance = np.inf

    for current_label in range(1, num_components + 1):
        y_pixels, x_pixels = np.where(labeled_mask == current_label)

        if x_pixels.size == 0:
            continue

        x_mean = np.mean(x_pixels)
        y_mean = np.mean(y_pixels)
        distance = np.hypot(x_mean - x_image_center, y_mean - y_image_center)

        if distance < best_distance:
            best_distance = distance
            best_label = current_label

    if best_label is None:
        raise RuntimeError("Could not identify the central connected component.")

    return labeled_mask == best_label


def estimate_component_center(component_mask: np.ndarray) -> tuple[float, float]:
    """
    Estimate the centroid of a connected component.

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.

    Returns
    -------
    x0 : float
        Centroid x coordinate in NumPy coordinates.
    y0 : float
        Centroid y coordinate in NumPy coordinates.
    """
    y_pixels, x_pixels = np.where(component_mask)

    if x_pixels.size == 0:
        raise RuntimeError("The selected component is empty.")

    x0 = np.mean(x_pixels)
    y0 = np.mean(y_pixels)

    return x0, y0


def pixel_is_selected(
    component_mask: np.ndarray,
    x: float,
    y: float,
) -> bool:
    """
    Test whether the nearest pixel to (x, y) belongs to the component mask.

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.
    x : float
        X coordinate in NumPy pixel coordinates.
    y : float
        Y coordinate in NumPy pixel coordinates.

    Returns
    -------
    bool
        True if the nearest pixel belongs to the selected component.
    """
    ny, nx = component_mask.shape

    ix = int(round(x))
    iy = int(round(y))

    if ix < 0 or ix >= nx or iy < 0 or iy >= ny:
        return False

    return bool(component_mask[iy, ix])


def ray_radius(
    component_mask: np.ndarray,
    x0: float,
    y0: float,
    theta_deg: float,
    step: float = 0.2,
    max_radius: float | None = None,
) -> float:
    """
    Measure the radius along one direction until the component mask ends.

    The angle is measured from the +x axis. Positive angles rotate
    counterclockwise on the displayed image.

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.
    x0 : float
        Center x coordinate in NumPy coordinates.
    y0 : float
        Center y coordinate in NumPy coordinates.
    theta_deg : float
        Direction angle in degrees.
    step : float, optional
        Radial step in pixels.
    max_radius : float, optional
        Maximum radius to test.

    Returns
    -------
    float
        Radius in pixels.
    """
    ny, nx = component_mask.shape

    if max_radius is None:
        max_radius = np.hypot(nx, ny)

    if not pixel_is_selected(component_mask, x0, y0):
        raise RuntimeError(
            "The estimated center is not inside the selected component. "
            "Try using --center median or inspect the mask."
        )

    theta_rad = np.deg2rad(theta_deg)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)

    last_good_radius = 0.0
    radius = step

    while radius <= max_radius:
        x = x0 + radius * cos_t
        y = y0 - radius * sin_t

        if not pixel_is_selected(component_mask, x, y):
            break

        last_good_radius = radius
        radius += step

    return last_good_radius


def symmetric_radii(
    component_mask: np.ndarray,
    x0: float,
    y0: float,
    theta_deg: float,
    radial_step: float = 0.2,
) -> tuple[float, float, float]:
    """
    Compute forward and opposite radii for one angle.

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.
    x0 : float
        Center x coordinate.
    y0 : float
        Center y coordinate.
    theta_deg : float
        Angle in degrees.
    radial_step : float, optional
        Radial step in pixels.

    Returns
    -------
    r_forward : float
        Radius at theta_deg.
    r_opposite : float
        Radius at theta_deg + 180 deg.
    score : float
        Symmetric score, defined as min(r_forward, r_opposite).
    """
    r_forward = ray_radius(
        component_mask=component_mask,
        x0=x0,
        y0=y0,
        theta_deg=theta_deg,
        step=radial_step,
    )

    r_opposite = ray_radius(
        component_mask=component_mask,
        x0=x0,
        y0=y0,
        theta_deg=theta_deg + 180.0,
        step=radial_step,
    )

    score = min(r_forward, r_opposite)

    return r_forward, r_opposite, score


def normalize_angle_180(theta_deg: float) -> float:
    """
    Normalize an angle to [0, 180).
    """
    return theta_deg % 180.0


def search_best_major_angle(
    component_mask: np.ndarray,
    x0: float,
    y0: float,
    angle_start: float = 0.0,
    angle_stop: float = 180.0,
    angle_step: float = 1.0,
    radial_step: float = 0.2,
) -> tuple[float, float, float, float]:
    """
    Search for the major-axis angle using a symmetric criterion.

    The selected angle maximizes min(r(theta), r(theta + 180)).
    Ties are broken using the larger sum r(theta) + r(theta + 180).

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.
    x0 : float
        Center x coordinate.
    y0 : float
        Center y coordinate.
    angle_start : float, optional
        Start angle in degrees.
    angle_stop : float, optional
        Stop angle in degrees.
    angle_step : float, optional
        Angular step in degrees.
    radial_step : float, optional
        Radial step in pixels.

    Returns
    -------
    best_theta : float
        Best angle in degrees.
    best_r1 : float
        Forward radius at the best angle.
    best_r2 : float
        Opposite radius at the best angle.
    best_score : float
        Symmetric score at the best angle.
    """
    if angle_step <= 0:
        raise ValueError("angle_step must be positive.")

    best_theta = None
    best_r1 = None
    best_r2 = None
    best_score = -np.inf
    best_sum = -np.inf

    angles = np.arange(angle_start, angle_stop, angle_step)

    for theta in angles:
        r1, r2, score = symmetric_radii(
            component_mask=component_mask,
            x0=x0,
            y0=y0,
            theta_deg=theta,
            radial_step=radial_step,
        )

        total = r1 + r2

        if (score > best_score) or (np.isclose(score, best_score) and total > best_sum):
            best_theta = normalize_angle_180(theta)
            best_r1 = r1
            best_r2 = r2
            best_score = score
            best_sum = total

    if best_theta is None:
        raise RuntimeError("Could not determine the major-axis angle.")

    return best_theta, best_r1, best_r2, best_score


def refine_major_angle(
    component_mask: np.ndarray,
    x0: float,
    y0: float,
    theta_coarse: float,
    coarse_step: float = 1.0,
    refine_window: float | None = None,
    refine_step: float = 0.1,
    radial_step: float = 0.2,
) -> tuple[float, float, float, float]:
    """
    Refine the major-axis angle around the coarse solution.

    Parameters
    ----------
    component_mask : numpy.ndarray
        Boolean mask of the selected component.
    x0 : float
        Center x coordinate.
    y0 : float
        Center y coordinate.
    theta_coarse : float
        Coarse major-axis angle.
    coarse_step : float, optional
        Coarse angular step.
    refine_window : float, optional
        Half-width of the refinement window.
    refine_step : float, optional
        Fine angular step.
    radial_step : float, optional
        Radial step in pixels.

    Returns
    -------
    best_theta : float
        Refined best angle in degrees.
    best_r1 : float
        Forward radius.
    best_r2 : float
        Opposite radius.
    best_score : float
        Symmetric score.
    """
    if refine_step <= 0:
        raise ValueError("refine_step must be positive.")

    if refine_window is None:
        refine_window = coarse_step

    angle_start = theta_coarse - refine_window
    angle_stop = theta_coarse + refine_window + 0.5 * refine_step

    return search_best_major_angle(
        component_mask=component_mask,
        x0=x0,
        y0=y0,
        angle_start=angle_start,
        angle_stop=angle_stop,
        angle_step=refine_step,
        radial_step=radial_step,
    )


def write_ds9_ellipse_region(
    region_file: str,
    x0_ds9: float,
    y0_ds9: float,
    semi_major: float,
    semi_minor: float,
    theta_ds9: float,
    color: str = "green",
) -> None:
    """
    Write a DS9 ellipse region file.

    Parameters
    ----------
    region_file : str
        Output region filename.
    x0_ds9 : float
        Center x coordinate in DS9 image coordinates.
    y0_ds9 : float
        Center y coordinate in DS9 image coordinates.
    semi_major : float
        Semi-major axis in pixels.
    semi_minor : float
        Semi-minor axis in pixels.
    theta_ds9 : float
        DS9 position angle in degrees.
    color : str, optional
        DS9 region color.
    """
    lines = [
        "# Region file format: DS9 version 4.1",
        (
            f"global color={color} dashlist=8 3 width=1 "
            'font="helvetica 10 normal roman" '
            "select=1 highlite=1 dash=0 fixed=0 "
            "edit=1 move=1 delete=1 include=1 source=1"
        ),
        "image",
        (
            f"ellipse({x0_ds9:.3f},{y0_ds9:.3f},"
            f"{semi_major:.3f},{semi_minor:.3f},{theta_ds9:.3f})"
        ),
    ]

    Path(region_file).write_text("\n".join(lines) + "\n")


def write_component_mask(
    output_file: str,
    component_mask: np.ndarray,
    header: fits.Header | None = None,
) -> None:
    """
    Write the selected central component as a FITS image.

    Parameters
    ----------
    output_file : str
        Output FITS filename.
    component_mask : numpy.ndarray
        Boolean component mask.
    header : astropy.io.fits.Header, optional
        Header copied from the input image.
    """
    hdu = fits.PrimaryHDU(data=component_mask.astype(np.uint8), header=header)
    hdu.header["IMTYPE"] = "component"
    hdu.writeto(output_file, overwrite=True)


def maincentralEllipse() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create a DS9 ellipse region for the central galaxy mask "
            "in a DESI maskbits FITS image."
        )
    )

    parser.add_argument(
        "maskbits",
        help="Input DESI maskbits FITS file.",
    )

    parser.add_argument(
        "region",
        help="Output DS9 ellipse region file.",
    )

    parser.add_argument(
        "--value",
        type=int,
        default=4096,
        help="Exact pixel value to use when --use-bit is not set. Default: 4096.",
    )

    parser.add_argument(
        "--use_bit",
        action="store_true",
        help="Use bit testing instead of exact pixel equality.",
    )

    parser.add_argument(
        "--bit",
        type=int,
        default=12,
        help="Bit number used when --use-bit is set. Default: 12.",
    )

    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Scale factor applied to the final ellipse axes. Default: 1.0.",
    )

    parser.add_argument(
        "--angle_step",
        type=float,
        default=1.0,
        help="Coarse angular step in degrees. Default: 1.0.",
    )

    parser.add_argument(
        "--refine",
        action="store_true",
        help="Enable local angular refinement.",
    )

    parser.add_argument(
        "--refine_window",
        type=float,
        default=None,
        help=(
            "Half-width of the refinement interval in degrees. "
            "If omitted, it is set equal to --angle-step."
        ),
    )

    parser.add_argument(
        "--refine_step",
        type=float,
        default=0.1,
        help="Fine angular step in degrees. Default: 0.1.",
    )

    parser.add_argument(
        "--radial_step",
        type=float,
        default=0.2,
        help="Radial step in pixels for ray tracing. Default: 0.2.",
    )

    parser.add_argument(
        "--component_out",
        default=None,
        help="Optional FITS file with the selected central component mask.",
    )

    parser.add_argument(
        "--color",
        default="green",
        help="DS9 region color. Default: green.",
    )

    args = parser.parse_args()

    central_ellipse(
        args.maskbits,
        args.region,
        args.value,
        args.use_bit,
        args.bit,
        args.scale,
        args.angle_step,
        args.refine,
        args.refine_window,
        args.refine_step,
        args.radial_step,
        args.component_out,
        args.color,
    )


def central_ellipse(
    maskbits,
    region,
    value=4096,
    use_bit=False,
    bit=12,
    scale=1.0,
    angle_step=1.0,
    refine=False,
    refine_window=None,
    refine_step=0.1,
    radial_step=0.2,
    component_out=None,
    color="green",
):

    if scale <= 0:
        raise ValueError("--scale must be positive.")

    image, header = read_first_image_hdu(maskbits)

    if image.ndim != 2:
        raise ValueError("The input image must be 2-dimensional.")

    selected_mask = build_selected_mask(
        image=image,
        target_value=value,
        use_bit=use_bit,
        bit=bit,
    )

    component_mask = select_central_component(selected_mask)

    x0_np, y0_np = estimate_component_center(component_mask)

    x0_ds9 = x0_np + 1.0
    y0_ds9 = y0_np + 1.0

    theta_major, r_major_1, r_major_2, score_major = search_best_major_angle(
        component_mask=component_mask,
        x0=x0_np,
        y0=y0_np,
        angle_start=0.0,
        angle_stop=180.0,
        angle_step=angle_step,
        radial_step=radial_step,
    )

    if refine:
        theta_major, r_major_1, r_major_2, score_major = refine_major_angle(
            component_mask=component_mask,
            x0=x0_np,
            y0=y0_np,
            theta_coarse=theta_major,
            coarse_step=angle_step,
            refine_window=refine_window,
            refine_step=refine_step,
            radial_step=radial_step,
        )

    semi_major = min(r_major_1, r_major_2)

    r_minor_1 = ray_radius(
        component_mask=component_mask,
        x0=x0_np,
        y0=y0_np,
        theta_deg=theta_major + 90.0,
        step=radial_step,
    )

    r_minor_2 = ray_radius(
        component_mask=component_mask,
        x0=x0_np,
        y0=y0_np,
        theta_deg=theta_major + 270.0,
        step=radial_step,
    )

    semi_minor = min(r_minor_1, r_minor_2)

    semi_major *= scale
    semi_minor *= scale

    # Convert ray-tracing angle to DS9 display angle.
    theta_ds9 = (-theta_major) % 180.0

    # Safety check: DS9 expects the first radius to be the major axis.
    if semi_minor > semi_major:
        semi_major, semi_minor = semi_minor, semi_major
        theta_ds9 = (theta_ds9 + 90.0) % 180.0

    write_ds9_ellipse_region(
        region_file=region,
        x0_ds9=x0_ds9,
        y0_ds9=y0_ds9,
        semi_major=semi_major,
        semi_minor=semi_minor,
        theta_ds9=theta_ds9,
        color=color,
    )

    if component_out is not None:
        write_component_mask(component_out, component_mask, header=header)

    print(f"Region written to: {region}")
    print(f"Center (NumPy)          = ({x0_np:.3f}, {y0_np:.3f})")
    print(f"Center (DS9)            = ({x0_ds9:.3f}, {y0_ds9:.3f})")
    print(f"Measured angle          = {theta_major:.3f} deg")
    print(f"DS9 angle               = {theta_ds9:.3f} deg")
    print(f"Major radius forward    = {r_major_1:.3f} pix")
    print(f"Major radius opposite   = {r_major_2:.3f} pix")
    print(f"Major symmetric score   = {score_major:.3f} pix")
    print(f"Minor radius forward    = {r_minor_1:.3f} pix")
    print(f"Minor radius opposite   = {r_minor_2:.3f} pix")
    print(f"Final semi-major axis   = {semi_major:.3f} pix")
    print(f"Final semi-minor axis   = {semi_minor:.3f} pix")

    if component_out is not None:
        print(f"Component mask written to: {component_out}")


if __name__ == "__main__":
    maincentralEllipse()
