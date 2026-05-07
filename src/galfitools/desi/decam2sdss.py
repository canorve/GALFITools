#!/usr/bin/env python3

# DECam (DESI South) ↔ SDSS transforms in g,r,z using linear color terms.
# Works with any available pair among (g,r), (r,z), (g,z).
# For each provided band B, it computes B_target = B_source + a + b * color,
# where 'color' is one of (g-r), (r-z), or (g-z) in the SOURCE system.

import numpy as np
import argparse

C_AA_S = 2.99792458e18  # speed of light [Å/s]
AB_ZP = 48.60  # AB zero-point [mag]

# ---------- SED (you may replace with a stellar library) ----------
def planck_flambda(wave_A, T_K):
    """Blackbody f_lambda (shape only), erg s^-1 cm^-2 Å^-1 (up to a constant)."""
    wave_A = np.asarray(wave_A, dtype=float)
    wave_cm = wave_A * 1e-8
    h = 6.62607015e-27  # erg s
    c = 2.99792458e10  # cm s^-1
    kB = 1.380649e-16  # erg K^-1
    x = (h * c) / (wave_cm * kB * T_K)
    Blam = (2.0 * h * c**2) / (wave_cm**5 * np.expm1(x))
    return Blam / 1e8  # per Å


# ---------- AB mag for photon-counting system ----------
def ab_mag_from_flambda(wave_A, trans, f_lambda):
    """Compute AB mag given wavelength [Å], throughput, and f_lambda."""
    wave = np.asarray(wave_A, dtype=float)
    S = np.asarray(trans, dtype=float)
    f = np.asarray(f_lambda, dtype=float)
    num = np.trapz(f * S * wave, wave)
    den = np.trapz(S * (C_AA_S / wave), wave)
    fnu = num / den
    return -2.5 * np.log10(fnu) - AB_ZP


# ---------- Load passbands (speclite) ----------
def load_filters_speclite():
    try:
        from speclite import filters
    except Exception as e:
        raise SystemExit(
            "speclite is required. Install with: pip install speclite"
        ) from e

    decam = filters.load_filters("decam2014-g", "decam2014-r", "decam2014-z")
    sdss = filters.load_filters("sdss2010-g", "sdss2010-r", "sdss2010-z")

    def as_angstrom(bp):
        w = bp.wavelength
        try:
            return w.to("angstrom").value  # Quantity -> Å
        except Exception:
            return np.asarray(w, dtype=float)  # already ndarray in Å

    F = {
        "decam_g": (as_angstrom(decam[0]), decam[0].response),
        "decam_r": (as_angstrom(decam[1]), decam[1].response),
        "decam_z": (as_angstrom(decam[2]), decam[2].response),
        "sdss_g": (as_angstrom(sdss[0]), sdss[0].response),
        "sdss_r": (as_angstrom(sdss[1]), sdss[1].response),
        "sdss_z": (as_angstrom(sdss[2]), sdss[2].response),
    }
    return F


# ---------- Linear fit with RMS ----------
def fit_linear(y, x):
    """Fit y = a + b x. Return (a, b, rms)."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    X = np.vstack([np.ones_like(x), x]).T
    a, b = np.linalg.lstsq(X, y, rcond=None)[0]
    rms = float(np.sqrt(np.mean((y - (a + b * x)) ** 2)))
    return float(a), float(b), rms


# ---------- Derive coefficients for all color bases ----------
def derive_coeffs(T_grid=(3500, 4000, 4500, 5000, 6000, 8000, 10000)):
    F = load_filters_speclite()
    wg, rg = F["decam_g"]
    wr, rr = F["decam_r"]
    wz, rz = F["decam_z"]
    wG, rG = F["sdss_g"]
    wR, rR = F["sdss_r"]
    wZ, rZ = F["sdss_z"]

    rows = []
    for T in T_grid:
        g_dec = ab_mag_from_flambda(wg, rg, planck_flambda(wg, T))
        r_dec = ab_mag_from_flambda(wr, rr, planck_flambda(wr, T))
        z_dec = ab_mag_from_flambda(wz, rz, planck_flambda(wz, T))
        g_sds = ab_mag_from_flambda(wG, rG, planck_flambda(wG, T))
        r_sds = ab_mag_from_flambda(wR, rR, planck_flambda(wR, T))
        z_sds = ab_mag_from_flambda(wZ, rZ, planck_flambda(wZ, T))
        rows.append((g_dec, r_dec, z_dec, g_sds, r_sds, z_sds))
    g_dec, r_dec, z_dec, g_sds, r_sds, z_sds = np.array(rows, dtype=float).T

    # Colors in each system
    c_gr_dec = g_dec - r_dec
    c_rz_dec = r_dec - z_dec
    c_gz_dec = g_dec - z_dec
    c_gr_sds = g_sds - r_sds
    c_rz_sds = r_sds - z_sds
    c_gz_sds = g_sds - z_sds

    def pack(delta, color):
        a, b, r = fit_linear(delta, color)
        return {"a": a, "b": b, "rms": r}

    coeffs = {
        "decam_to_sdss": {
            "gr": {
                "g": pack(g_sds - g_dec, c_gr_dec),
                "r": pack(r_sds - r_dec, c_gr_dec),
                "z": pack(z_sds - z_dec, c_gr_dec),
            },
            "rz": {
                "g": pack(g_sds - g_dec, c_rz_dec),
                "r": pack(r_sds - r_dec, c_rz_dec),
                "z": pack(z_sds - z_dec, c_rz_dec),
            },
            "gz": {
                "g": pack(g_sds - g_dec, c_gz_dec),
                "r": pack(r_sds - r_dec, c_gz_dec),
                "z": pack(z_sds - z_dec, c_gz_dec),
            },
        },
        "sdss_to_decam": {
            "gr": {
                "g": pack(g_dec - g_sds, c_gr_sds),
                "r": pack(r_dec - r_sds, c_gr_sds),
                "z": pack(z_dec - z_sds, c_gr_sds),
            },
            "rz": {
                "g": pack(g_dec - g_sds, c_rz_sds),
                "r": pack(r_dec - r_sds, c_rz_sds),
                "z": pack(z_dec - z_sds, c_rz_sds),
            },
            "gz": {
                "g": pack(g_dec - g_sds, c_gz_sds),
                "r": pack(r_dec - r_sds, c_gz_sds),
                "z": pack(z_dec - z_sds, c_gz_sds),
            },
        },
        "T_grid": tuple(T_grid),
    }
    return coeffs


# ---------- Apply transform (safe with missing bands) ----------
def apply_transform(direction, mags_in, coeffs, color_base=None):
    """
    direction: 'decam_to_sdss' or 'sdss_to_decam'
    mags_in: dict with provided source mags, e.g., {'g':16.99, 'z':15.20}
    color_base: 'gr','rz','gz' or None to auto-choose from the provided pair.
    Returns (out_dict, base_used, rms_dict). Only bands present in mags_in are output.
    """
    have = {k for k, v in mags_in.items() if v is not None and np.isfinite(v)}

    # choose base
    bases = []
    if {"g", "r"} <= have:
        bases.append("gr")
    if {"r", "z"} <= have:
        bases.append("rz")
    if {"g", "z"} <= have:
        bases.append("gz")

    if color_base is None:
        if not bases:
            raise ValueError("Provide at least one pair: (g,r), (r,z), or (g,z).")
        base = bases[0]
    else:
        if color_base not in ("gr", "rz", "gz"):
            raise ValueError("color_base must be one of: gr, rz, gz.")
        if not set(color_base) <= have:
            raise ValueError(f"color_base='{color_base}' requires {tuple(color_base)}.")
        base = color_base

    c = float(mags_in[base[0]] - mags_in[base[1]])

    out = {}
    rms = {}
    for band in ("g", "r", "z"):
        src = mags_in.get(band, None)
        if src is None or not np.isfinite(src):
            continue  # only transform bands you actually provided
        a = coeffs[direction][base][band]["a"]
        b = coeffs[direction][base][band]["b"]
        out[band] = float(src) + a + b * c
        rms[band] = coeffs[direction][base][band]["rms"]

    return out, base, rms


# ---------- CLI ----------
def maindecam2sdss():
    ap = argparse.ArgumentParser(
        description="DECam ↔ SDSS transforms in g,r,z using any color pair."
    )
    ap.add_argument(
        "--temps",
        type=float,
        nargs="+",
        default=[3500, 4000, 4500, 5000, 6000, 8000, 10000],
        help="Blackbody temperatures [K] to fit color terms.",
    )
    ap.add_argument(
        "--print-coeffs",
        action="store_true",
        help="Print all coefficients and RMS, then exit.",
    )
    ap.add_argument(
        "--color-base",
        choices=["gr", "rz", "gz"],
        default=None,
        help="Force a specific color base (default: auto from provided pair).",
    )

    sub = ap.add_subparsers(dest="cmd")

    a = sub.add_parser("decam2sdss", help="Convert DECam → SDSS.")
    a.add_argument("--g", type=float)
    a.add_argument("--r", type=float)
    a.add_argument("--z", type=float)

    b = sub.add_parser("sdss2decam", help="Convert SDSS → DECam.")
    b.add_argument("--g", type=float)
    b.add_argument("--r", type=float)
    b.add_argument("--z", type=float)

    args = ap.parse_args()
    coeffs = derive_coeffs(tuple(args.temps))

    if args.print_coeffs or args.cmd is None:
        print("# Coefficients (Δ = target - source) and RMS [mag]:")
        for direction in ("decam_to_sdss", "sdss_to_decam"):
            print(f"[{direction}]")
            for base in ("gr", "rz", "gz"):
                parts = []
                for band in ("g", "r", "z"):
                    c = coeffs[direction][base][band]
                    parts.append(
                        f"{band}: {c['a']:+.5f} + {c['b']:+.5f} * {base} (rms={c['rms']:.4f})"
                    )
                print("  " + base + ": " + " | ".join(parts))
        print(f"T grid [K]: {coeffs['T_grid']}")
        if args.cmd is None:
            return

    if args.cmd == "decam2sdss":
        mags_in = {"g": args.g, "r": args.r, "z": args.z}
        out, base, rms = apply_transform(
            "decam_to_sdss", mags_in, coeffs, args.color_base
        )
        print(f"# Used color base: {base}")
        for band in ("g", "r", "z"):
            if band in out:
                print(
                    f"{band}_SDSS = {out[band]:.4f}   (fit rms ≈ {rms[band]:.4f} mag)"
                )

    elif args.cmd == "sdss2decam":
        mags_in = {"g": args.g, "r": args.r, "z": args.z}
        out, base, rms = apply_transform(
            "sdss_to_decam", mags_in, coeffs, args.color_base
        )
        print(f"# Used color base: {base}")
        for band in ("g", "r", "z"):
            if band in out:
                print(
                    f"{band}_DECam = {out[band]:.4f}   (fit rms ≈ {rms[band]:.4f} mag)"
                )


if __name__ == "__main__":
    maindecam2sdss()
