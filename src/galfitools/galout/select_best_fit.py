#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


FLOAT_TEXT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
NUMBER_PATTERN = re.compile(FLOAT_TEXT)
WRAPPED_NUMBER_PATTERN = re.compile(
    rf"(?P<opening>[\[{{*]?)(?P<number>{FLOAT_TEXT})(?P<closing>[\]}}*]?)"
)
MODEL_PATTERN = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*(.*)$")
CHI_NU_PATTERN = re.compile(rf"^\s*Chi\^2/nu\s*=\s*({FLOAT_TEXT})\s*$")
CHI2_PATTERN = re.compile(
    rf"^\s*Chi\^2\s*=\s*({FLOAT_TEXT})\s*,\s*ndof\s*=\s*(\d+)\s*$"
)
INPUT_IMAGE_PATTERN = re.compile(r"^\s*Input image\s*:\s*(.+?)\s*$")
OUTPUT_IMAGE_PATTERN = re.compile(r"^\s*Output image\s*:\s*(.+?)\s*$")
INIT_FILE_PATTERN = re.compile(r"^\s*Init\. par\. file\s*:\s*(.+?)\s*$")
RESTART_FILE_PATTERN = re.compile(r"^\s*Restart file\s*:\s*(.+?)\s*$")
IMAGE_REGION_PATTERN = re.compile(
    r"^(?P<image>.+)\["
    r"(?P<xmin>\d+)\s*:\s*(?P<xmax>\d+)\s*,\s*"
    r"(?P<ymin>\d+)\s*:\s*(?P<ymax>\d+)"
    r"\]\s*$"
)
SEPARATOR_PATTERN = re.compile(r"^-{20,}\s*$")
DEFAULT_SEPARATOR = "-" * 77


@dataclass(frozen=True)
class LogHeader:
    """Information extracted from the selected fit.log entry."""

    input_image: str
    output_image: str
    fit_region: tuple[int, int, int, int]
    init_file: str | None
    restart_file: str | None
    chi2: float | None
    ndof: int | None
    chi_nu: float


@dataclass(frozen=True)
class Parameter:
    """One selected GALFIT parameter and its fitting state."""

    index: int
    name: str
    value: float
    error: float
    toggle: int


@dataclass(frozen=True)
class Component:
    """One GALFIT model component from the selected entry."""

    number: int
    model: str
    parameters: tuple[Parameter, ...]


@dataclass(frozen=True)
class ParameterLine:
    """Mapping from fit.log parameter order to GALFIT input syntax."""

    galfit_number: str
    parameter_indices: tuple[int, ...]
    description: str


MODEL_PARAMETER_NAMES: dict[str, tuple[str, ...]] = {
    "sersic": (
        "x",
        "y",
        "magnitude",
        "r_e",
        "n",
        "axis_ratio",
        "position_angle",
    ),
    "expdisk": (
        "x",
        "y",
        "magnitude",
        "r_s",
        "axis_ratio",
        "position_angle",
    ),
    "devauc": (
        "x",
        "y",
        "magnitude",
        "r_e",
        "axis_ratio",
        "position_angle",
    ),
    "gaussian": (
        "x",
        "y",
        "magnitude",
        "fwhm",
        "axis_ratio",
        "position_angle",
    ),
    "moffat": (
        "x",
        "y",
        "magnitude",
        "fwhm",
        "power_law",
        "axis_ratio",
        "position_angle",
    ),
    "ferrer": (
        "x",
        "y",
        "magnitude",
        "r_out",
        "alpha",
        "beta",
        "axis_ratio",
        "position_angle",
    ),
    "psf": (
        "x",
        "y",
        "magnitude",
    ),
    "sky": (
        "sky",
        "dsky_dx",
        "dsky_dy",
    ),
}


MODEL_OUTPUT_LINES: dict[str, tuple[ParameterLine, ...]] = {
    "sersic": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "R_e (effective radius) [pix]"),
        ParameterLine("5", (5,), "Sersic index n"),
        ParameterLine("9", (6,), "Axis ratio (b/a)"),
        ParameterLine("10", (7,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "expdisk": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "Scale length [pix]"),
        ParameterLine("9", (5,), "Axis ratio (b/a)"),
        ParameterLine("10", (6,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "devauc": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "R_e (effective radius) [pix]"),
        ParameterLine("9", (5,), "Axis ratio (b/a)"),
        ParameterLine("10", (6,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "gaussian": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "FWHM [pix]"),
        ParameterLine("9", (5,), "Axis ratio (b/a)"),
        ParameterLine("10", (6,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "moffat": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "FWHM [pix]"),
        ParameterLine("5", (5,), "Power-law index"),
        ParameterLine("9", (6,), "Axis ratio (b/a)"),
        ParameterLine("10", (7,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "ferrer": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
        ParameterLine("4", (4,), "Outer truncation radius [pix]"),
        ParameterLine("5", (5,), "Alpha"),
        ParameterLine("6", (6,), "Beta"),
        ParameterLine("9", (7,), "Axis ratio (b/a)"),
        ParameterLine("10", (8,), "Position angle (PA) [deg: Up=0, Left=90]"),
    ),
    "psf": (
        ParameterLine("1", (1, 2), "Position x, y"),
        ParameterLine("3", (3,), "Integrated magnitude"),
    ),
    "sky": (
        ParameterLine("1", (1,), "Sky background at center [ADUs]"),
        ParameterLine("2", (2,), "dsky/dx [ADUs/pix]"),
        ParameterLine("3", (3,), "dsky/dy [ADUs/pix]"),
    ),
}


MODEL_EXTRA_LINES: dict[str, tuple[str, ...]] = {
    "sersic": (
        " 6) 0.0000      0          # -----",
        " 7) 0.0000      0          # -----",
        " 8) 0.0000      0          # -----",
    ),
}


def split_entries(log_text: str) -> tuple[str, list[str]]:
    """
    Split fit.log into entries.

    Entry lengths and numbers of model components may vary. Separator lines
    are removed here and reconstructed when the selected log is written.
    """
    entries: list[str] = []
    current_entry: list[str] = []
    separator = DEFAULT_SEPARATOR

    for line in log_text.splitlines():
        if SEPARATOR_PATTERN.fullmatch(line.strip()):
            separator = line.strip()
            continue

        if line.lstrip().startswith("Input image"):
            if current_entry:
                entries.append("\n".join(current_entry).strip())
            current_entry = [line]
        elif current_entry:
            current_entry.append(line)

    if current_entry:
        entries.append("\n".join(current_entry).strip())

    return separator, entries


def extract_chi_nu(entry: str, entry_number: int) -> float:
    """Extract the displayed reduced chi-square value from one entry."""
    for line in entry.splitlines():
        match = CHI_NU_PATTERN.match(line)
        if match:
            return float(match.group(1))

    raise ValueError(f"Could not find 'Chi^2/nu' in entry {entry_number}.")


def select_lowest_chi_nu(entries: Sequence[str]) -> tuple[int, str, float]:
    """
    Select the first entry with the minimum displayed Chi^2/nu value.

    A strict '<' comparison intentionally preserves the first entry when two
    or more entries have exactly the same displayed reduced chi-square.
    """
    best_index = -1
    best_entry = ""
    best_chi_nu = float("inf")

    for entry_number, entry in enumerate(entries, start=1):
        chi_nu = extract_chi_nu(entry, entry_number)

        if chi_nu < best_chi_nu:
            best_index = entry_number
            best_entry = entry
            best_chi_nu = chi_nu

    if best_index < 0:
        raise ValueError("No valid fit.log entries were found.")

    return best_index, best_entry, best_chi_nu


def extract_number_tokens(text: str) -> list[tuple[float, int]]:
    """Extract values and convert fit.log wrappers into GALFIT toggles."""
    tokens: list[tuple[float, int]] = []

    for match in WRAPPED_NUMBER_PATTERN.finditer(text):
        opening = match.group("opening")
        closing = match.group("closing")

        if "*" in opening or "*" in closing:
            raise ValueError(f"Starred parameter found in line: {text}")

        if opening == "[" or closing == "]":
            toggle = 0
        elif opening == "{" or closing == "}":
            toggle = 2
        else:
            toggle = 1

        tokens.append((float(match.group("number")), toggle))

    return tokens


def extract_numbers(text: str) -> list[float]:
    """Extract floating-point values while ignoring wrappers."""
    return [float(number) for number in NUMBER_PATTERN.findall(text)]


def get_parameter_names(model: str, count: int) -> tuple[str, ...]:
    """Return the expected names for one supported GALFIT component."""
    names = MODEL_PARAMETER_NAMES.get(model)

    if names is None:
        raise ValueError(f"Model '{model}' is not supported for GALFIT-file output.")

    if len(names) != count:
        raise ValueError(
            f"Model '{model}' has {count} values in fit.log, but the script "
            f"expects {len(names)}."
        )

    return names


def parse_components(entry: str) -> list[Component]:
    """Parse all model components and errors from the selected entry."""
    lines = entry.splitlines()
    components: list[Component] = []
    line_index = 0

    while line_index < len(lines):
        model_match = MODEL_PATTERN.match(lines[line_index])

        if model_match is None:
            line_index += 1
            continue

        model = model_match.group(1).lower()
        fitted_tokens = extract_number_tokens(model_match.group(2))

        error_line_index = line_index + 1
        while error_line_index < len(lines) and not lines[error_line_index].strip():
            error_line_index += 1

        if error_line_index >= len(lines):
            raise ValueError(f"Missing uncertainty line below model '{model}'.")

        errors = extract_numbers(lines[error_line_index])

        if not errors:
            raise ValueError(
                f"Could not read the uncertainty line for model '{model}'."
            )

        if len(fitted_tokens) < len(errors):
            raise ValueError(
                f"Model '{model}' has {len(fitted_tokens)} fitted values but "
                f"{len(errors)} uncertainties."
            )

        # The sky line contains two reference coordinates before its three
        # actual fitted parameters. Keep the final N values, where N is the
        # number of uncertainties on the line below.
        fitted_tokens = fitted_tokens[-len(errors) :]
        names = get_parameter_names(model, len(fitted_tokens))

        parameters = tuple(
            Parameter(
                index=parameter_index,
                name=name,
                value=value,
                error=abs(error),
                toggle=toggle,
            )
            for parameter_index, (name, (value, toggle), error) in enumerate(
                zip(names, fitted_tokens, errors),
                start=1,
            )
        )

        components.append(
            Component(
                number=len(components) + 1,
                model=model,
                parameters=parameters,
            )
        )

        line_index = error_line_index + 1

    if not components:
        raise ValueError("No GALFIT components were found in the selected entry.")

    return components


def extract_header(entry: str, chi_nu: float) -> LogHeader:
    """Extract GALFIT control values available in the selected entry."""
    input_value: str | None = None
    output_image: str | None = None
    init_file: str | None = None
    restart_file: str | None = None
    chi2: float | None = None
    ndof: int | None = None

    for line in entry.splitlines():
        if match := INPUT_IMAGE_PATTERN.match(line):
            input_value = match.group(1).strip()
        elif match := OUTPUT_IMAGE_PATTERN.match(line):
            output_image = match.group(1).strip()
        elif match := INIT_FILE_PATTERN.match(line):
            init_file = match.group(1).strip()
        elif match := RESTART_FILE_PATTERN.match(line):
            restart_file = match.group(1).strip()
        elif match := CHI2_PATTERN.match(line):
            chi2 = float(match.group(1))
            ndof = int(match.group(2))

    if input_value is None:
        raise ValueError("Could not find 'Input image' in the selected entry.")

    if output_image is None:
        raise ValueError("Could not find 'Output image' in the selected entry.")

    region_match = IMAGE_REGION_PATTERN.match(input_value)
    if region_match is None:
        raise ValueError(
            "Could not extract the fitting region from 'Input image'. "
            "Use --fit-region XMIN XMAX YMIN YMAX."
        )

    return LogHeader(
        input_image=region_match.group("image").strip(),
        output_image=output_image,
        fit_region=(
            int(region_match.group("xmin")),
            int(region_match.group("xmax")),
            int(region_match.group("ymin")),
            int(region_match.group("ymax")),
        ),
        init_file=init_file,
        restart_file=restart_file,
        chi2=chi2,
        ndof=ndof,
        chi_nu=chi_nu,
    )


def write_selected_log(
    output_file: Path,
    entry: str,
    separator: str,
) -> None:
    """Write only the selected fit.log entry."""
    text = f"{separator}\n\n{entry.strip()}\n\n{separator}\n"
    output_file.write_text(text, encoding="utf-8")


def format_parameter_value(value: float, parameter_name: str) -> str:
    """Format a selected value for a GALFIT input file."""
    if parameter_name in {"sky", "dsky_dx", "dsky_dy"}:
        return f"{value:.8e}"
    return f"{value:.6f}"


def write_component(component: Component) -> list[str]:
    """Convert one selected component to a GALFIT input block."""
    line_definitions = MODEL_OUTPUT_LINES.get(component.model)

    if line_definitions is None:
        raise ValueError(
            f"No GALFIT output mapping is defined for model '{component.model}'."
        )

    by_index = {parameter.index: parameter for parameter in component.parameters}
    lines = [
        f"# Component number: {component.number}",
        f" 0) {component.model:<22s} # Component type",
    ]

    for line_definition in line_definitions:
        if line_definition.galfit_number == "9":
            lines.extend(MODEL_EXTRA_LINES.get(component.model, ()))

        selected = [by_index[index] for index in line_definition.parameter_indices]
        values = " ".join(
            format_parameter_value(parameter.value, parameter.name)
            for parameter in selected
        )
        toggles = " ".join(str(parameter.toggle) for parameter in selected)
        errors = ", ".join(f"{parameter.error:.6g}" for parameter in selected)

        lines.append(
            f"{line_definition.galfit_number:>2s}) "
            f"{values} {toggles}  # {line_definition.description}; "
            f"GALFIT error(s): {errors}"
        )

    lines.append(" Z) 0                      # Skip this model? (yes=1, no=0)")
    return lines


def write_galfit_file(
    output_file: Path,
    components: Sequence[Component],
    header: LogHeader,
    source_log: Path,
    selected_entry_number: int,
    sigma_image: str,
    psf_image: str,
    psf_sampling: int,
    mask_image: str,
    constraints_file: str,
    convolution_box: tuple[int, int],
    zeropoint: float,
    plate_scale: tuple[float, float],
    display_type: str,
    run_option: int,
) -> None:
    """Write a complete GALFIT file using the selected fitted values."""
    xmin, xmax, ymin, ymax = header.fit_region
    conv_x, conv_y = convolution_box
    scale_x, scale_y = plate_scale

    lines = [
        f"# GALFIT file generated from: {source_log}",
        f"# Selected fit.log entry: {selected_entry_number}",
        f"# Input parameter file: {header.init_file or 'unknown'}",
        f"# Restart file: {header.restart_file or 'unknown'}",
        f"# Chi^2/nu = {header.chi_nu:.12g}",
    ]

    if header.chi2 is not None and header.ndof is not None:
        lines.append(f"# Chi^2 = {header.chi2:.12g}, ndof = {header.ndof}")

    lines.extend(
        [
            "",
            "================================================================================",
            "# IMAGE and GALFIT CONTROL PARAMETERS",
            f"A) {header.input_image:<45s} # Input data image (FITS file)",
            f"B) {header.output_image:<45s} # Output data image block",
            f"C) {sigma_image:<45s} # Sigma image",
            f"D) {psf_image:<45s} # Input PSF image",
            f"E) {psf_sampling:<45d} # PSF fine sampling factor",
            f"F) {mask_image:<45s} # Bad pixel mask",
            f"G) {constraints_file:<45s} # Parameter constraints file",
            f"H) {xmin} {xmax} {ymin} {ymax:<31d} # Image region to fit",
            f"I) {conv_x} {conv_y:<39d} # Convolution box size (x y)",
            f"J) {zeropoint:<45.6f} # Magnitude photometric zeropoint",
            f"K) {scale_x:.6f} {scale_y:<36.6f} # Plate scale [arcsec/pixel]",
            f"O) {display_type:<45s} # Display type",
            f"P) {run_option:<45d} # 0=optimize, 1=model, 2=imgblock, 3=subcomps",
            "",
            "# INITIAL FITTING PARAMETERS",
            "# Values and fitting toggles are copied from the selected fit.",
            "# [] in fit.log -> 0; {} -> 2; ordinary value -> 1.",
            "",
        ]
    )

    for component in components:
        lines.extend(write_component(component))
        lines.append("")

    lines.append(
        "================================================================================"
    )
    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select the first GALFIT fit.log entry with the lowest displayed "
            "Chi^2/nu, then write a one-entry fit.log and a GALFIT input file."
        )
    )

    parser.add_argument("fit_log", type=Path, help="Input GALFIT fit.log file")
    parser.add_argument(
        "-o",
        "--log-output",
        type=Path,
        default=Path("best_fit.log"),
        help="Selected one-entry log (default: best_fit.log)",
    )
    parser.add_argument(
        "--galfit-output",
        type=Path,
        default=Path("best_fit.init"),
        help="Selected GALFIT parameter file (default: best_fit.init)",
    )

    parser.add_argument("--input-image", help="Override GALFIT A)")
    parser.add_argument("--output-image", help="Override GALFIT B)")
    parser.add_argument(
        "--fit-region",
        nargs=4,
        type=int,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Override GALFIT H)",
    )

    parser.add_argument(
        "--sigma-image",
        default="sigma.fits",
        help="GALFIT C) value (default: sigma.fits)",
    )
    parser.add_argument(
        "--psf-image",
        default="psf.fits",
        help="GALFIT D) value (default: psf.fits)",
    )
    parser.add_argument(
        "--psf-sampling",
        type=int,
        default=1,
        help="GALFIT E) value (default: 1)",
    )
    parser.add_argument(
        "--mask-image",
        default="mask.fts",
        help="GALFIT F) value (default: mask.fts)",
    )
    parser.add_argument(
        "--constraints-file",
        default="constraints.txt",
        help="GALFIT G) value (default: constraints.txt)",
    )
    parser.add_argument(
        "--convolution-box",
        nargs=2,
        type=int,
        default=(100, 100),
        metavar=("X", "Y"),
        help="GALFIT I) values (default: 100 100)",
    )
    parser.add_argument(
        "--zeropoint",
        type=float,
        default=22.5,
        help="GALFIT J) value (default: 22.5)",
    )
    parser.add_argument(
        "--plate-scale",
        nargs=2,
        type=float,
        default=(0.262, 0.262),
        metavar=("DX", "DY"),
        help="GALFIT K) values (default: 0.262 0.262)",
    )
    parser.add_argument(
        "--display-type",
        choices=("regular", "curses", "both"),
        default="regular",
        help="GALFIT O) value (default: regular)",
    )
    parser.add_argument(
        "--run-option",
        type=int,
        choices=(0, 1, 2, 3),
        default=2,
        help="GALFIT P) value (default: 2)",
    )

    return parser


def mainSelectBestFitlog() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not args.fit_log.is_file():
        parser.error(f"Input file does not exist: {args.fit_log}")

    try:
        log_text = args.fit_log.read_text(encoding="utf-8")
        separator, entries = split_entries(log_text)

        if not entries:
            raise ValueError(f"No GALFIT entries were found in {args.fit_log}.")

        selected_number, selected_entry, selected_chi_nu = select_lowest_chi_nu(entries)
        header = extract_header(selected_entry, selected_chi_nu)
        components = parse_components(selected_entry)

        if args.fit_region is not None:
            header = LogHeader(
                input_image=header.input_image,
                output_image=header.output_image,
                fit_region=tuple(args.fit_region),
                init_file=header.init_file,
                restart_file=header.restart_file,
                chi2=header.chi2,
                ndof=header.ndof,
                chi_nu=header.chi_nu,
            )

        header = LogHeader(
            input_image=args.input_image or header.input_image,
            output_image=args.output_image or header.output_image,
            fit_region=header.fit_region,
            init_file=header.init_file,
            restart_file=header.restart_file,
            chi2=header.chi2,
            ndof=header.ndof,
            chi_nu=header.chi_nu,
        )

        write_selected_log(args.log_output, selected_entry, separator)
        write_galfit_file(
            output_file=args.galfit_output,
            components=components,
            header=header,
            source_log=args.fit_log,
            selected_entry_number=selected_number,
            sigma_image=args.sigma_image,
            psf_image=args.psf_image,
            psf_sampling=args.psf_sampling,
            mask_image=args.mask_image,
            constraints_file=args.constraints_file,
            convolution_box=tuple(args.convolution_box),
            zeropoint=args.zeropoint,
            plate_scale=tuple(args.plate_scale),
            display_type=args.display_type,
            run_option=args.run_option,
        )

    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"Entries read       : {len(entries)}")
    print(f"Selected entry     : {selected_number}")
    print(f"Chi^2/nu           : {selected_chi_nu:.12g}")
    print(f"Init. par. file    : {header.init_file or 'unknown'}")
    if header.chi2 is not None:
        print(f"Chi^2              : {header.chi2:.12g}")
    print(f"Selected log       : {args.log_output}")
    print(f"GALFIT parameter   : {args.galfit_output}")


if __name__ == "__main__":
    mainSelectBestFitlog()
