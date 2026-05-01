#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
import time
from io import BytesIO
from typing import Iterable, List, Tuple

import numpy as np
import requests
from astropy.io import fits
from requests.adapters import HTTPAdapter
from requests.exceptions import (
    ChunkedEncodingError,
    ConnectionError,
    HTTPError,
    Timeout,
)
from urllib3.util.retry import Retry


CUTOUT_URL = "https://www.legacysurvey.org/viewer/cutout.fits"
COADD_PSF_URL = "https://www.legacysurvey.org/viewer/coadd-psf/"

# Errors that are likely to be temporary during a large download.
RETRYABLE_EXCEPTIONS = (ChunkedEncodingError, ConnectionError, Timeout)


def read_radec_csv(path: str) -> List[Tuple[str, float, float]]:
    """
    Read a CSV file with RA and DEC columns.

    The CSV file must contain columns named ``ra`` and ``dec``. These names are
    case-insensitive. If an ``objid`` column exists, it is used as the object
    identifier. Otherwise, sequential identifiers are assigned.

    Parameters
    ----------
    path : str
        Input CSV filename.

    Returns
    -------
    list of tuple
        Each tuple contains ``(objid, ra, dec)``.
    """
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(path)

    with p.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header.")

        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        if "ra" not in field_map or "dec" not in field_map:
            raise ValueError("CSV must contain columns named 'ra' and 'dec'.")

        objid_key = field_map.get("objid")
        ra_key = field_map["ra"]
        dec_key = field_map["dec"]

        out: List[Tuple[str, float, float]] = []
        for count, row in enumerate(reader, start=1):
            objid = str(row[objid_key]).strip() if objid_key else f"obj{count:04d}"
            out.append((objid, float(row[ra_key]), float(row[dec_key])))

    return out


def configure_session() -> requests.Session:
    """Create a requests session with retries for common HTTP status errors."""
    session = requests.Session()

    retry = Retry(
        total=3,
        connect=3,
        read=3,
        status=3,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    return session


def fetch_fits(
    session: requests.Session,
    url: str,
    params: dict,
    timeout: int = 120,
    retries: int = 5,
    retry_wait: float = 3.0,
    label: str = "FITS file",
) -> fits.HDUList:
    """
    Download a FITS file with retries.

    This function catches incomplete reads and temporary connection failures.
    These errors can occur when a remote server closes a large download before
    all bytes have arrived.

    Parameters
    ----------
    session : requests.Session
        Active HTTP session.
    url : str
        URL of the FITS service.
    params : dict
        Query parameters for the request.
    timeout : int, optional
        Request timeout in seconds.
    retries : int, optional
        Maximum number of attempts.
    retry_wait : float, optional
        Initial waiting time between attempts, in seconds. The wait time grows
        linearly with the attempt number.
    label : str, optional
        Short label printed in the error messages.

    Returns
    -------
    astropy.io.fits.HDUList
        Downloaded FITS file.
    """
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            response = session.get(url, params=params, timeout=timeout)
            response.raise_for_status()

            # response.content can raise ChunkedEncodingError if the server sends
            # fewer bytes than announced in the response header.
            data = response.content
            return fits.open(BytesIO(data), memmap=False)

        except HTTPError as err:
            last_error = err
            status = err.response.status_code if err.response is not None else None

            # HTTP 4xx errors are usually not fixed by retrying, except 429.
            if status is not None and 400 <= status < 500 and status != 429:
                raise RuntimeError(
                    f"{label}: HTTP error {status}. Parameters: {params}"
                ) from err

        except RETRYABLE_EXCEPTIONS as err:
            last_error = err

        except OSError as err:
            # astropy may raise OSError if the downloaded bytes are not a valid
            # or complete FITS file.
            last_error = err

        if attempt < retries:
            wait = retry_wait * attempt
            print(
                f"Warning: {label} download failed on attempt "
                f"{attempt}/{retries}: {last_error}. Retrying in {wait:.1f} s."
            )
            time.sleep(wait)

    raise RuntimeError(
        f"{label}: failed after {retries} attempts. Parameters: {params}. "
        f"Last error: {last_error}"
    ) from last_error


def write_primary(data, header, outpath: pathlib.Path) -> None:
    """Write data and header as a primary FITS HDU."""
    outpath.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(outpath, overwrite=True)


def find_first_image_hdu(hdul: fits.HDUList) -> fits.ImageHDU | fits.PrimaryHDU:
    """Return the first HDU that contains image data."""
    for hdu in hdul:
        if hdu.data is not None:
            return hdu

    raise RuntimeError("Could not find an image HDU in the returned FITS.")


def find_hdu_by_name_or_index(
    hdul: fits.HDUList,
    names: set[str],
    fallback_index: int | None = None,
) -> fits.ImageHDU | fits.PrimaryHDU | None:
    """Return an HDU by extension name, with an optional index fallback."""
    for hdu in hdul:
        hdu_name = hdu.name.strip().upper()
        if hdu.data is not None and hdu_name in names:
            return hdu

    if fallback_index is not None and len(hdul) > fallback_index:
        hdu = hdul[fallback_index]
        if hdu.data is not None:
            return hdu

    return None


def split_cutout_products(
    hdul: fits.HDUList,
    require_maskbits: bool = True,
) -> tuple[fits.PrimaryHDU, fits.PrimaryHDU, fits.PrimaryHDU | None]:
    """Return the image, inverse-variance, and optional maskbits HDUs.

    The Legacy Survey FITS cutout service may return several extensions.
    The usual order is image, inverse variance, and maskbits when both
    ``invvar`` and ``maskbits`` are requested. This function first searches by
    HDU name and then falls back to the expected extension order.
    """
    img_hdu = find_first_image_hdu(hdul)
    inv_hdu = find_hdu_by_name_or_index(
        hdul,
        names={"INVVAR", "IVAR", "WEIGHT"},
        fallback_index=1,
    )
    maskbits_hdu = find_hdu_by_name_or_index(
        hdul,
        names={"MASKBITS", "MASK", "BITMASK"},
        fallback_index=2,
    )

    if inv_hdu is None:
        raise RuntimeError("Could not find an invvar HDU in the returned FITS.")

    if require_maskbits and maskbits_hdu is None:
        raise RuntimeError("Could not find a maskbits HDU in the returned FITS.")

    img = fits.PrimaryHDU(data=img_hdu.data, header=img_hdu.header)
    inv = fits.PrimaryHDU(data=inv_hdu.data, header=inv_hdu.header)

    maskbits = None
    if maskbits_hdu is not None:
        maskbits = fits.PrimaryHDU(data=maskbits_hdu.data, header=maskbits_hdu.header)
        maskbits.header["IMTYPE"] = "maskbits"

    return img, inv, maskbits


def convert_to_sigma(image_file: pathlib.Path) -> None:
    """
    Convert an inverse-variance DESI image into a sigma image for GALFIT.

    DESI inverse variance is assumed to be ``invvar = 1 / sigma**2``.
    Therefore, the sigma image is computed as ``sigma = 1 / sqrt(invvar)``.
    Invalid or non-positive inverse-variance values are assigned a large sigma.
    """
    with fits.open(image_file) as hdu:
        invvar = np.asarray(hdu[0].data).astype(np.float64, copy=False)

        sigma = np.full(invvar.shape, np.nan, dtype=np.float64)
        valid = np.isfinite(invvar) & (invvar > 0.0)
        sigma[valid] = 1.0 / np.sqrt(invvar[valid])
        sigma[~valid] = 1e6

        hdu_new = fits.PrimaryHDU(sigma, hdu[0].header)
        hdu_new.header["IMTYPE"] = "sigma"
        hdu_new.header["BUNIT"] = "nanomaggy"

        output_file = image_file.with_name(image_file.name.replace("invvar", "sigma"))
        hdu_new.writeto(output_file, overwrite=True)


def write_failure_log(failures: list[dict[str, str]], outdir: pathlib.Path) -> None:
    """Write a CSV file with the downloads that failed."""
    if not failures:
        return

    log_path = outdir / "download_failures.csv"
    with log_path.open("w", newline="") as f:
        fieldnames = ["objid", "ra", "dec", "band", "stage", "error"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(failures)

    print(f"Failed downloads were written to: {log_path}")


def main_downloadDesi() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Download DESI Legacy Survey image, invvar, maskbits, and PSF cutouts. "
            "The invvar image is also converted to a GALFIT sigma image."
        )
    )
    ap.add_argument("csv", help="Input CSV with columns ra and dec. Optional: objid.")
    ap.add_argument("--outdir", default="legacy_cutouts", help="Output directory")
    ap.add_argument(
        "--layer", default="ls-dr10", help="Viewer layer, e.g. ls-dr10 or ls-dr9"
    )
    ap.add_argument("--size", type=int, default=256, help="Cutout size in pixels")
    ap.add_argument(
        "--pixscale", type=float, default=0.262, help="Arcsec/pixel for cutouts"
    )
    ap.add_argument("--bands", default="grz", help="Bands to download, e.g. grz")
    ap.add_argument(
        "--subimage",
        action="store_true",
        help=(
            "Add the subimage flag. This uses the fixed brick grid and includes "
            "the inverse-variance image."
        ),
    )
    ap.add_argument(
        "--no-maskbits",
        action="store_true",
        help="Do not request or write the DESI maskbits cutout.",
    )
    ap.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Timeout in seconds for each request.",
    )
    ap.add_argument(
        "--retries",
        type=int,
        default=5,
        help="Maximum number of attempts for each FITS download.",
    )
    ap.add_argument(
        "--retry-wait",
        type=float,
        default=3.0,
        help="Initial waiting time in seconds between retries.",
    )
    ap.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop the program when a download fails after all retries.",
    )

    args = ap.parse_args()

    targets = read_radec_csv(args.csv)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bands: Iterable[str] = list(args.bands.strip())
    failures: list[dict[str, str]] = []

    with configure_session() as session:
        for i, (objid, ra, dec) in enumerate(targets, start=1):
            obj_dir = outdir / objid
            base_prefix = f"{objid}_ra{ra:.6f}_dec{dec:.6f}"

            print(f"Processing {i}/{len(targets)}: {objid}  RA={ra:.6f}  DEC={dec:.6f}")

            for band in bands:
                band_dir = obj_dir / band
                band_dir.mkdir(parents=True, exist_ok=True)

                cut_params = {
                    "ra": ra,
                    "dec": dec,
                    "layer": args.layer,
                    "size": args.size,
                    "pixscale": args.pixscale,
                    "bands": band,
                    "invvar": 1,
                }
                if not args.no_maskbits:
                    cut_params["maskbits"] = 1
                if args.subimage:
                    cut_params["subimage"] = 1

                try:
                    label = f"cutout objid={objid} band={band}"
                    with fetch_fits(
                        session,
                        CUTOUT_URL,
                        cut_params,
                        timeout=args.timeout,
                        retries=args.retries,
                        retry_wait=args.retry_wait,
                        label=label,
                    ) as hdul:
                        img_hdu, inv_hdu, maskbits_hdu = split_cutout_products(
                            hdul, require_maskbits=not args.no_maskbits
                        )

                    img_path = band_dir / f"{base_prefix}_{band}_img.fits"
                    inv_path = band_dir / "invvar.fits"
                    maskbits_path = band_dir / "maskbits.fits"

                    img_hdu.writeto(img_path, overwrite=True)
                    inv_hdu.writeto(inv_path, overwrite=True)
                    if maskbits_hdu is not None:
                        maskbits_hdu.writeto(maskbits_path, overwrite=True)
                    convert_to_sigma(inv_path)

                except Exception as err:
                    message = str(err)
                    print(
                        f"Error: cutout failed for objid={objid}, band={band}: {message}"
                    )
                    failures.append(
                        {
                            "objid": objid,
                            "ra": f"{ra:.8f}",
                            "dec": f"{dec:.8f}",
                            "band": band,
                            "stage": "cutout",
                            "error": message,
                        }
                    )
                    if args.stop_on_error:
                        write_failure_log(failures, outdir)
                        raise
                    continue

                psf_params = {"ra": ra, "dec": dec, "layer": args.layer, "bands": band}

                try:
                    label = f"PSF objid={objid} band={band}"
                    with fetch_fits(
                        session,
                        COADD_PSF_URL,
                        psf_params,
                        timeout=args.timeout,
                        retries=args.retries,
                        retry_wait=args.retry_wait,
                        label=label,
                    ) as psf_hdul:
                        psf_data_hdu = (
                            psf_hdul[0] if psf_hdul[0].data is not None else psf_hdul[1]
                        )
                        psf_path = band_dir / "psf.fits"
                        write_primary(psf_data_hdu.data, psf_data_hdu.header, psf_path)

                except Exception as err:
                    message = str(err)
                    print(
                        f"Error: PSF failed for objid={objid}, band={band}: {message}"
                    )
                    failures.append(
                        {
                            "objid": objid,
                            "ra": f"{ra:.8f}",
                            "dec": f"{dec:.8f}",
                            "band": band,
                            "stage": "psf",
                            "error": message,
                        }
                    )
                    if args.stop_on_error:
                        write_failure_log(failures, outdir)
                        raise
                    continue

    write_failure_log(failures, outdir)

    if failures:
        print(f"Download finished with {len(failures)} failed item(s).")
        return 1

    print("Download done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main_downloadDesi())
