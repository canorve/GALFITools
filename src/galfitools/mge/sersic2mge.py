#! /usr/bin/env python3


from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List
import os

import numpy as np

try:
    from scipy.optimize import nnls

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


from galfitools.galin.galfit import (
    GalComps,
    Galfit,
    GalHead,
    GetRadAng,
    SelectGal,
    conver2Sersic,
    conver2Edge,
    numComps,
    galPrintComp,
    galPrintHeader,
    galPrintSky,
)


@dataclass
class SersicParameters:
    """
    Parameters of a Sérsic component.

    Parameters
    ----------
    x0 : float
        X position in pixels.
    y0 : float
        Y position in pixels.
    magnitude : float
        Integrated magnitude.
    re_pix : float
        Effective radius in pixels.
    n : float
        Sérsic index.
    axis_ratio : float
        Axis ratio q = b/a.
    pa_deg : float
        Position angle in degrees using the GALFIT convention.
    zeropoint : float
        Photometric zeropoint.
    """

    x0: float
    y0: float
    magnitude: float
    re_pix: float
    n: float
    axis_ratio: float
    pa_deg: float
    zeropoint: float


@dataclass
class GaussianFluxComponent:
    """
    Gaussian component represented by total flux.

    Parameters
    ----------
    flux : float
        Total flux of the Gaussian.
    sigma_pix : float
        Sigma in pixels.
    axis_ratio : float
        Axis ratio q = b/a.
    """

    flux: float
    sigma_pix: float
    axis_ratio: float


@dataclass
class GalfitGaussianComponent:
    """
    Gaussian component in GALFIT format.

    Parameters
    ----------
    x0 : float
        X position in pixels.
    y0 : float
        Y position in pixels.
    magnitude : float
        Integrated magnitude.
    fwhm_pix : float
        FWHM in pixels.
    axis_ratio : float
        Axis ratio q = b/a.
    pa_deg : float
        Position angle in degrees.
    """

    x0: float
    y0: float
    magnitude: float
    fwhm_pix: float
    axis_ratio: float
    pa_deg: float


def sersic_bn(n: float) -> float:
    """
    Compute the Sérsic b_n coefficient.

    Parameters
    ----------
    n : float
        Sérsic index.

    Returns
    -------
    float
        The b_n coefficient.
    """
    if n <= 0:
        raise ValueError("Sérsic index n must be positive.")

    if np.isclose(n, 0.5, atol=1e-12):
        return math.log(2.0)

    return (
        2.0 * n
        - 1.0 / 3.0
        + 4.0 / (405.0 * n)
        + 46.0 / (25515.0 * n**2)
        + 131.0 / (1148175.0 * n**3)
        - 2194697.0 / (30690717750.0 * n**4)
    )


def total_flux_from_magnitude(magnitude: float, zeropoint: float) -> float:
    """
    Convert integrated magnitude into total flux.

    Parameters
    ----------
    magnitude : float
        Integrated magnitude.
    zeropoint : float
        Photometric zeropoint.

    Returns
    -------
    float
        Total flux.
    """
    return 10.0 ** (-0.4 * (magnitude - zeropoint))


def magnitude_from_flux(flux: float, zeropoint: float) -> float:
    """
    Convert total flux into magnitude.

    Parameters
    ----------
    flux : float
        Total flux.
    zeropoint : float
        Photometric zeropoint.

    Returns
    -------
    float
        Integrated magnitude.
    """
    if flux <= 0:
        raise ValueError("Flux must be positive.")
    return zeropoint - 2.5 * math.log10(flux)


def sersic_ie_from_total_flux(params: SersicParameters) -> float:
    """
    Compute I_e from the total flux of an elliptical Sérsic model.

    Parameters
    ----------
    params : SersicParameters
        Sérsic parameters.

    Returns
    -------
    float
        Intensity at the effective radius.
    """
    if params.re_pix <= 0:
        raise ValueError("re_pix must be positive.")
    if params.axis_ratio <= 0 or params.axis_ratio > 1:
        raise ValueError("axis_ratio must satisfy 0 < q <= 1.")

    total_flux = total_flux_from_magnitude(params.magnitude, params.zeropoint)
    b_n = sersic_bn(params.n)

    denominator = (
        2.0
        * math.pi
        * params.axis_ratio
        * params.n
        * math.exp(b_n)
        * params.re_pix**2
        * math.gamma(2.0 * params.n)
        / (b_n ** (2.0 * params.n))
    )
    return total_flux / denominator


def sersic_profile(r: np.ndarray, params: SersicParameters) -> np.ndarray:
    """
    Evaluate the Sérsic profile.

    Parameters
    ----------
    r : np.ndarray
        Radius array in pixels.
    params : SersicParameters
        Sérsic parameters.

    Returns
    -------
    np.ndarray
        Intensity profile.
    """
    ie = sersic_ie_from_total_flux(params)
    b_n = sersic_bn(params.n)
    x = np.maximum(r / params.re_pix, 0.0)
    return ie * np.exp(-b_n * (x ** (1.0 / params.n) - 1.0))


def gaussian_unit_flux_profile(
    r: np.ndarray,
    sigma_pix: float,
    axis_ratio: float,
) -> np.ndarray:
    """
    Evaluate a Gaussian profile normalized to unit total flux.

    Parameters
    ----------
    r : np.ndarray
        Radius array in pixels.
    sigma_pix : float
        Gaussian sigma in pixels.
    axis_ratio : float
        Axis ratio q = b/a.

    Returns
    -------
    np.ndarray
        Gaussian intensity profile with total flux equal to 1.
    """
    norm = 1.0 / (2.0 * math.pi * axis_ratio * sigma_pix**2)
    return norm * np.exp(-0.5 * (r / sigma_pix) ** 2)


def choose_sigma_grid(
    re_pix: float,
    n_gaussians: int,
    min_sigma_factor: float = 0.02,
    max_sigma_factor: float = 8.0,
) -> np.ndarray:
    """
    Create a logarithmic sigma grid.

    Parameters
    ----------
    re_pix : float
        Effective radius in pixels.
    n_gaussians : int
        Number of Gaussian components.
    min_sigma_factor : float, optional
        Minimum sigma relative to R_e.
    max_sigma_factor : float, optional
        Maximum sigma relative to R_e.

    Returns
    -------
    np.ndarray
        Sigma grid in pixels.
    """
    if re_pix <= 0:
        raise ValueError("re_pix must be positive.")
    if n_gaussians < 1:
        raise ValueError("n_gaussians must be at least 1.")
    if min_sigma_factor <= 0 or max_sigma_factor <= 0:
        raise ValueError("Sigma factors must be positive.")
    if min_sigma_factor >= max_sigma_factor:
        raise ValueError("min_sigma_factor must be smaller than max_sigma_factor.")

    return np.geomspace(
        min_sigma_factor * re_pix, max_sigma_factor * re_pix, n_gaussians
    )


def fit_sersic_with_flux_conservation(
    params: SersicParameters,
    n_gaussians: int = 6,
    r_max_factor: float = 10.0,
    n_samples: int = 800,
    min_sigma_factor: float = 0.02,
    max_sigma_factor: float = 8.0,
) -> List[GaussianFluxComponent]:
    """
    Approximate a Sérsic profile with Gaussian components conserving total flux.

    The fit is done in terms of total Gaussian fluxes F_j.
    After the fit, the solution is renormalized so that the sum of all
    Gaussian fluxes matches exactly the total Sérsic flux.

    Parameters
    ----------
    params : SersicParameters
        Sérsic parameters.
    n_gaussians : int, optional
        Number of Gaussian components.
    r_max_factor : float, optional
        Maximum radius in units of R_e.
    n_samples : int, optional
        Number of radial samples.
    min_sigma_factor : float, optional
        Minimum sigma relative to R_e.
    max_sigma_factor : float, optional
        Maximum sigma relative to R_e.

    Returns
    -------
    list of GaussianFluxComponent
        Gaussian components with flux-conserving normalization.
    """
    total_flux = total_flux_from_magnitude(params.magnitude, params.zeropoint)

    if np.isclose(params.n, 0.5, atol=1e-12):
        sigma_pix = params.re_pix / math.sqrt(2.0 * math.log(2.0))
        return [
            GaussianFluxComponent(
                flux=total_flux,
                sigma_pix=sigma_pix,
                axis_ratio=params.axis_ratio,
            )
        ]

    r = np.geomspace(1e-4, r_max_factor * params.re_pix, n_samples)
    target = sersic_profile(r, params)

    sigma_grid = choose_sigma_grid(
        re_pix=params.re_pix,
        n_gaussians=n_gaussians,
        min_sigma_factor=min_sigma_factor,
        max_sigma_factor=max_sigma_factor,
    )

    design_matrix = np.column_stack(
        [
            gaussian_unit_flux_profile(r, sigma, params.axis_ratio)
            for sigma in sigma_grid
        ]
    )

    weights = 1.0 / np.sqrt(r)
    weighted_matrix = design_matrix * weights[:, None]
    weighted_target = target * weights

    if SCIPY_AVAILABLE:
        fluxes, _ = nnls(weighted_matrix, weighted_target)
    else:
        fluxes, _, _, _ = np.linalg.lstsq(weighted_matrix, weighted_target, rcond=None)
        fluxes = np.clip(fluxes, 0.0, None)

    flux_sum = float(np.sum(fluxes))
    if flux_sum <= 0:
        raise RuntimeError("The Gaussian fit produced zero total flux.")

    # Exact flux conservation
    fluxes *= total_flux / flux_sum

    components = [
        GaussianFluxComponent(
            flux=float(flux),
            sigma_pix=float(sigma),
            axis_ratio=params.axis_ratio,
        )
        for flux, sigma in zip(fluxes, sigma_grid)
        if flux > 0.0
    ]

    return components


def sigma_to_fwhm(sigma_pix: float) -> float:
    """
    Convert sigma to FWHM.

    Parameters
    ----------
    sigma_pix : float
        Sigma in pixels.

    Returns
    -------
    float
        FWHM in pixels.
    """
    return 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma_pix


def gaussian_profile_from_flux(
    r: np.ndarray,
    component: GaussianFluxComponent,
) -> np.ndarray:
    """
    Evaluate one Gaussian profile from total flux.

    Parameters
    ----------
    r : np.ndarray
        Radius array in pixels.
    component : GaussianFluxComponent
        Gaussian component.

    Returns
    -------
    np.ndarray
        Intensity profile.
    """
    return component.flux * gaussian_unit_flux_profile(
        r,
        component.sigma_pix,
        component.axis_ratio,
    )


def combined_gaussian_profile(
    r: np.ndarray,
    components: List[GaussianFluxComponent],
) -> np.ndarray:
    """
    Evaluate the summed Gaussian profile.

    Parameters
    ----------
    r : np.ndarray
        Radius array in pixels.
    components : list of GaussianFluxComponent
        Gaussian components.

    Returns
    -------
    np.ndarray
        Summed intensity profile.
    """
    model = np.zeros_like(r, dtype=float)
    for component in components:
        model += gaussian_profile_from_flux(r, component)
    return model


def convert_to_galfit_gaussians(
    params: SersicParameters,
    components: List[GaussianFluxComponent],
) -> List[GalfitGaussianComponent]:
    """
    Convert Gaussian flux components to GALFIT Gaussian parameters.

    Parameters
    ----------
    params : SersicParameters
        Original Sérsic parameters.
    components : list of GaussianFluxComponent
        Gaussian components.

    Returns
    -------
    list of GalfitGaussianComponent
        Components in GALFIT format.
    """
    galfit_components: List[GalfitGaussianComponent] = []

    for component in components:
        mag = magnitude_from_flux(component.flux, params.zeropoint)
        fwhm = sigma_to_fwhm(component.sigma_pix)

        galfit_components.append(
            GalfitGaussianComponent(
                x0=params.x0,
                y0=params.y0,
                magnitude=mag,
                fwhm_pix=fwhm,
                axis_ratio=component.axis_ratio,
                pa_deg=params.pa_deg,
            )
        )

    return galfit_components


def format_galfit_gaussian_block(
    components: List[GalfitGaussianComponent],
    start_index: int = 1,
    fit_position: bool = False,
    fit_magnitude: bool = True,
    fit_fwhm: bool = True,
    fit_axis_ratio: bool = False,
    fit_pa: bool = False,
) -> str:
    """
    Format Gaussian components as GALFIT text blocks.

    Parameters
    ----------
    components : list of GalfitGaussianComponent
        Components in GALFIT format.
    start_index : int, optional
        Starting component number.
    fit_position : bool, optional
        Whether x and y are free.
    fit_magnitude : bool, optional
        Whether magnitude is free.
    fit_fwhm : bool, optional
        Whether FWHM is free.
    fit_axis_ratio : bool, optional
        Whether axis ratio is free.
    fit_pa : bool, optional
        Whether PA is free.

    Returns
    -------
    str
        GALFIT block text.
    """
    pos_flag = "1 1" if fit_position else "0 0"
    mag_flag = "1" if fit_magnitude else "0"
    fwhm_flag = "1" if fit_fwhm else "0"
    q_flag = "1" if fit_axis_ratio else "0"
    pa_flag = "1" if fit_pa else "0"

    lines: List[str] = []

    for i, comp in enumerate(components, start=start_index):
        lines.extend(
            [
                f"# Component number: {i}",
                " 0) gaussian                 #  Component type",
                f" 1) {comp.x0:.4f} {comp.y0:.4f} {pos_flag}  #  Position x, y",
                f" 3) {comp.magnitude:.6f}     {mag_flag}      #  Integrated magnitude",
                f" 4) {comp.fwhm_pix:.6f}      {fwhm_flag}      #  FWHM [pix]",
                f" 9) {comp.axis_ratio:.6f}      {q_flag}      #  Axis ratio (b/a)",
                f"10) {comp.pa_deg:.6f}      {pa_flag}      #  Position angle (PA) [deg: Up=0, Left=90]",
                " Z) 0                      #  Skip this model in output image?  (yes=1, no=0)",
                "",
            ]
        )

    return "\n".join(lines)


def total_magnitude_of_gaussians(
    components: List[GaussianFluxComponent],
    zeropoint: float,
) -> float:
    """
    Compute the total integrated magnitude of all Gaussian components.

    Parameters
    ----------
    components : list of GaussianFluxComponent
        Gaussian components.
    zeropoint : float
        Photometric zeropoint.

    Returns
    -------
    float
        Total magnitude.
    """
    total_flux = sum(component.flux for component in components)
    return magnitude_from_flux(total_flux, zeropoint)


def print_summary(
    params: SersicParameters,
    components: List[GaussianFluxComponent],
) -> None:
    """
    Print a summary of the flux consistency.

    Parameters
    ----------
    params : SersicParameters
        Original Sérsic parameters.
    components : list of GaussianFluxComponent
        Gaussian components.
    """
    original_flux = total_flux_from_magnitude(params.magnitude, params.zeropoint)
    gaussian_flux = sum(component.flux for component in components)

    original_mag = params.magnitude
    gaussian_mag = magnitude_from_flux(gaussian_flux, params.zeropoint)

    print("Flux consistency check")
    print("----------------------")
    print(f"Original Sérsic magnitude : {original_mag:.6f}")
    print(f"Gaussian total magnitude  : {gaussian_mag:.6f}")
    print(f"Original total flux       : {original_flux:.6e}")
    print(f"Gaussian total flux       : {gaussian_flux:.6e}")
    print(f"Flux ratio G/S            : {gaussian_flux / original_flux:.12f}")
    print()


def Sersic2mge(args) -> None:

    galfit = Galfit(args.GalfitFile)

    galhead = galfit.ReadHead()
    galcomps = galfit.ReadComps()
    galsky = galfit.ReadSky()

    num_comp = 1  # assumes it will only convert the first galaxy
    dis = 3

    # convert all exp, gaussian and de vaucouleurs to Sersic format
    comps = conver2Sersic(galcomps)

    comps.Flux = 10 ** ((-comps.Mag) / 2.5)

    comps = SelectGal(comps, dis, num_comp)

    # printing output file
    fout = open(args.output, "w")

    filename = galhead.outimage
    root, extension = os.path.splitext(filename)
    newname = root + "-mge.fits"
    galhead.outimage = newname

    galPrintHeader(fout, galhead)

    index = 0

    for index, item in enumerate(comps.N):
        conver2mge(
            comps,
            args.numgauss,
            args.rmax,
            args.nsamples,
            args.minsigma,
            args.maxsigma,
            galhead.mgzpt,
            index,
            fout,
        )

    galPrintSky(fout, index + 1, galsky)


def conver2mge(
    galcomps, numgauss, rmax, nsamples, minsigma, maxsigma, zp, index, fout
) -> None:
    """Documentar


    Parameters
    ----------


    """
    params = SersicParameters(
        x0=galcomps.PosX[index],
        y0=galcomps.PosY[index],
        magnitude=galcomps.Mag[index],
        re_pix=galcomps.Rad[index],
        n=galcomps.Exp[index],
        axis_ratio=galcomps.AxRat[index],
        pa_deg=galcomps.PosAng[index],
        zeropoint=zp,
    )

    n_gaussians = 6

    components = fit_sersic_with_flux_conservation(
        params=params,
        n_gaussians=numgauss,
        r_max_factor=rmax,
        n_samples=nsamples,
        min_sigma_factor=minsigma,
        max_sigma_factor=maxsigma,
    )

    print_summary(params, components)

    galfit_components = convert_to_galfit_gaussians(params, components)
    galfit_text = format_galfit_gaussian_block(
        galfit_components,
        start_index=1,
        fit_position=True,
        fit_magnitude=True,
        fit_fwhm=True,
        fit_axis_ratio=True,
        fit_pa=True,
    )

    print(galfit_text)

    fout.write(galfit_text)
