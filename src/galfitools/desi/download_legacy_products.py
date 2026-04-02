#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import pathlib
from io import BytesIO
from typing import Iterable, List, Tuple

import requests
from astropy.io import fits
import numpy as np
import sys
import os


CUTOUT_URL = "https://www.legacysurvey.org/viewer/cutout.fits"
COADD_PSF_URL = "https://www.legacysurvey.org/viewer/coadd-psf/"


def read_radec_csv(path: str) -> List[Tuple[float, float]]:
    """
    Reads a CSV that contains at least ra, dec columns (case-insensitive).
    Extra columns are ignored.
    """
    p = pathlib.Path(path)
    if not p.exists():
        raise FileNotFoundError(path)

    count = 1
    with p.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("CSV has no header.")

        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        if "ra" not in field_map or "dec" not in field_map:
            raise ValueError(
                "CSV must contain columns named 'ra' and 'dec' (any case)."
            )

        if "objid" in field_map:
            objid_key = field_map["objid"]

        ra_key = field_map["ra"]
        dec_key = field_map["dec"]

        out: List[Tuple[float, float]] = []
        for row in reader:

            if "objid" in field_map:
                out.append(
                    (str(row[objid_key]), float(row[ra_key]), float(row[dec_key]))
                )
            else:
                out.append((str(row[count]), float(row[ra_key]), float(row[dec_key])))
                count = count + 1
        return out


def fetch_fits(
    session: requests.Session, url: str, params: dict, timeout: int = 120
) -> fits.HDUList:
    r = session.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return fits.open(BytesIO(r.content))


def write_primary(data, header, outpath: pathlib.Path) -> None:
    outpath.parent.mkdir(parents=True, exist_ok=True)
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(outpath, overwrite=True)


def split_cutout_invvar(hdul: fits.HDUList) -> tuple[fits.PrimaryHDU, fits.PrimaryHDU]:
    """Return (image, invvar) as PrimaryHDUs."""
    img_hdu = hdul[0]

    inv_hdu = None
    if len(hdul) >= 2:
        inv_hdu = hdul[1]
    else:
        for h in hdul:
            if h.name.strip().upper() in {"INVVAR", "IVAR", "WEIGHT"}:
                inv_hdu = h
                break

    if inv_hdu is None:
        raise RuntimeError("Could not find an invvar HDU in the returned FITS.")

    img = fits.PrimaryHDU(data=img_hdu.data, header=img_hdu.header)
    inv = fits.PrimaryHDU(data=inv_hdu.data, header=inv_hdu.header)
    return img, inv


def main_downloadDesi() -> int:
    ap = argparse.ArgumentParser(
        description="Download image, invvar, mask images from DESI and converts invvar to sigma image "
    )
    ap.add_argument(
        "csv", help="Input CSV with at least columns: ra, dec. optional: objid"
    )
    ap.add_argument("--outdir", default="legacy_cutouts", help="Output directory")
    ap.add_argument(
        "--layer", default="ls-dr10", help="Viewer layer, e.g. ls-dr10 or ls-dr9"
    )
    ap.add_argument(
        "--size", type=int, default=256, help="Cutout size in pixels (square)"
    )
    ap.add_argument(
        "--pixscale", type=float, default=0.262, help="Arcsec/pixel for cutouts"
    )
    ap.add_argument("--bands", default="grz", help="Bands to download, e.g. grz")
    ap.add_argument(
        "--subimage",
        action="store_true",
        help="If set, adds 'subimage' flag (no resampling; fixed brick grid; includes invvar).",
    )
    args = ap.parse_args()

    targets = read_radec_csv(args.csv)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bands: Iterable[str] = list(args.bands.strip())

    with requests.Session() as session:
        for i, (objid, ra, dec) in enumerate(targets, start=1):
            # Folder name: obj0001, obj0002, ...
            # obj_folder = f"obj{i:04d}"
            # obj_folder = f"obj{i}"
            obj_folder = objid
            obj_dir = outdir / obj_folder

            # File prefix: obj1_ra..._dec..._<band>_<type>.fits
            obj_id = f"obj{i}"
            # base_prefix = f"{obj_id}_ra{ra:.6f}_dec{dec:.6f}"
            base_prefix = f"{objid}_ra{ra:.6f}_dec{dec:.6f}"

            for b in bands:
                band_dir = obj_dir / b
                band_dir.mkdir(parents=True, exist_ok=True)

                # 1) image + invvar cutout
                cut_params = {
                    "ra": ra,
                    "dec": dec,
                    "layer": args.layer,
                    "size": args.size,
                    "pixscale": args.pixscale,
                    "bands": b,
                    "invvar": 1,
                }
                if args.subimage:
                    cut_params["subimage"] = 1

                hdul = fetch_fits(session, CUTOUT_URL, cut_params)
                img_hdu, inv_hdu = split_cutout_invvar(hdul)

                img_path = band_dir / f"{base_prefix}_{b}_img.fits"
                # inv_path = band_dir / f"{base_prefix}_{b}_invvar.fits"
                inv_path = band_dir / f"invvar.fits"
                img_hdu.writeto(img_path, overwrite=True)
                inv_hdu.writeto(inv_path, overwrite=True)

                # 2) coadd PSF
                psf_params = {"ra": ra, "dec": dec, "layer": args.layer, "bands": b}
                psf_hdul = fetch_fits(session, COADD_PSF_URL, psf_params)
                psf_data_hdu = (
                    psf_hdul[0] if psf_hdul[0].data is not None else psf_hdul[1]
                )

                # psf_path = band_dir / f"{base_prefix}_{b}_psf.fits"
                psf_path = band_dir / f"psf.fits"
                write_primary(psf_data_hdu.data, psf_data_hdu.header, psf_path)

                # converting to sigma image:
                convert_to_sigma(inv_path)

    print("Download done.")

    return 0


def convert_to_sigma(image_file):
    """
    Convert an inverse variance DESI image (invvar = 1/sigma^2)
    to a sigma image for GALFIT (sigma = 1/sqrt(invvar)).
    """
    with fits.open(image_file) as hdu:
        # NumPy 2.0-safe: asarray allows a copy when unavoidable
        invvar = np.asarray(hdu[0].data).astype(np.float64, copy=False)

        sigma = np.full(invvar.shape, np.nan, dtype=np.float64)

        m = np.isfinite(invvar) & (invvar > 0.0)
        sigma[m] = 1.0 / np.sqrt(invvar[m])

        # If GALFIT does not tolerate NaNs, replace invalid pixels by a large sigma:
        sigma[~m] = 1e6

        hdu_new = fits.PrimaryHDU(sigma, hdu[0].header)
        hdu_new.header["IMTYPE"] = "sigma"
        hdu_new.header["BUNIT"] = "nanomaggy"

        output_file = image_file.with_name(image_file.name.replace("invvar", "sigma"))

        hdu_new.writeto(output_file, overwrite=True)


# if __name__ == "__main__":
#    raise SystemExit(main())
