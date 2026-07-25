#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


FLOAT_TEXT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
NUMBER_PATTERN = re.compile(FLOAT_TEXT)
WRAPPED_NUMBER_PATTERN = re.compile(
    rf"(?P<opening>[\[{{*]?)(?P<number>{FLOAT_TEXT})(?P<closing>[\]}}*]?)"
)
MODEL_PATTERN = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*(.*)$")
STARRED_PARAMETER_PATTERN = re.compile(r"\*[^*\n]+\*")
INPUT_IMAGE_PATTERN = re.compile(r"^\s*Input image\s*:\s*(.+?)\s*$")
OUTPUT_IMAGE_PATTERN = re.compile(r"^\s*Output image\s*:\s*(.+?)\s*$")
IMAGE_REGION_PATTERN = re.compile(
    r"^(?P<image>.+)\["
    r"(?P<xmin>\d+)\s*:\s*(?P<xmax>\d+)\s*,\s*"
    r"(?P<ymin>\d+)\s*:\s*(?P<ymax>\d+)"
    r"\]\s*$"
)


@dataclass(frozen=True)
class ParameterKey:
    """Unique identification of one GALFIT model parameter."""

    component: int
    model: str
    index: int
    name: str


@dataclass(frozen=True)
class Measurement:
    """One fitted value, GALFIT uncertainty, and fit toggle."""

    value: float
    error: float
    toggle: int


@dataclass(frozen=True)
class ParameterStatistics:
    """Statistics calculated from all GALFIT entries."""

    key: ParameterKey
    toggle: int
    number_of_fits: int
    mean: float
    formal_error_on_mean: float
    scatter_standard_deviation: float
    scatter_error_on_mean: float
    combined_error: float
    minimum: float
    maximum: float


@dataclass(frozen=True)
class LogHeader:
    """Header information extracted from fit.log."""

    input_image: str
    output_image: str
    fit_region: tuple[int, int, int, int]


@dataclass(frozen=True)
class ParameterLine:
    """Definition of a line in a GALFIT component block."""

    galfit_number: str
    parameter_indices: tuple[int, ...]
    description: str


# The indices refer to the order in which parameters appear in fit.log.
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


# Lines inserted because they do not appear in the compact fit.log output.
MODEL_EXTRA_LINES: dict[str, tuple[str, ...]] = {
    "sersic": (
        " 6) 0.0000      0          # -----",
        " 7) 0.0000      0          # -----",
        " 8) 0.0000      0          # -----",
    ),
}


def extract_number_tokens(text: str) -> list[tuple[float, int]]:
    """Extract values and infer GALFIT toggles from [], {}, or no wrapper."""
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
    """Extract floating-point values while ignoring brackets and braces."""
    return [float(number) for number in NUMBER_PATTERN.findall(text)]


def split_entries(log_text: str) -> list[str]:
    """Split fit.log into entries beginning with an Input image line."""
    entries: list[str] = []
    current_entry: list[str] = []

    for line in log_text.splitlines():
        if line.lstrip().startswith("Input image"):
            if current_entry:
                entries.append("\n".join(current_entry))
            current_entry = [line]
        elif current_entry:
            current_entry.append(line)

    if current_entry:
        entries.append("\n".join(current_entry))

    return entries


def extract_log_header(entry: str) -> LogHeader:
    """Extract A, B, and H values from one fit.log entry."""
    input_value: str | None = None
    output_image: str | None = None

    for line in entry.splitlines():
        input_match = INPUT_IMAGE_PATTERN.match(line)
        if input_match:
            input_value = input_match.group(1).strip()

        output_match = OUTPUT_IMAGE_PATTERN.match(line)
        if output_match:
            output_image = output_match.group(1).strip()

    if input_value is None:
        raise ValueError("Could not find 'Input image' in a fit.log entry.")

    if output_image is None:
        raise ValueError("Could not find 'Output image' in a fit.log entry.")

    region_match = IMAGE_REGION_PATTERN.match(input_value)
    if region_match is None:
        raise ValueError(
            "Could not extract the fitting region from the Input image line. "
            "Use --fit-region XMIN XMAX YMIN YMAX."
        )

    input_image = region_match.group("image").strip()
    fit_region = (
        int(region_match.group("xmin")),
        int(region_match.group("xmax")),
        int(region_match.group("ymin")),
        int(region_match.group("ymax")),
    )

    return LogHeader(
        input_image=input_image,
        output_image=output_image,
        fit_region=fit_region,
    )


def get_parameter_names(model: str, number_of_parameters: int) -> tuple[str, ...]:
    """Return parameter names and reject unsupported output structures."""
    names = MODEL_PARAMETER_NAMES.get(model)

    if names is None:
        raise ValueError(
            f"Model '{model}' is not supported for GALFIT-file output. "
            "Add its parameter names and GALFIT line mapping to the script."
        )

    if len(names) != number_of_parameters:
        raise ValueError(
            f"Model '{model}' has {number_of_parameters} parameters in fit.log, "
            f"but the script expects {len(names)}."
        )

    return names


def parse_entry(
    entry: str,
    entry_number: int,
) -> dict[ParameterKey, Measurement]:
    """Extract each component value, uncertainty, and fit toggle."""
    if STARRED_PARAMETER_PATTERN.search(entry):
        raise ValueError(
            f"Entry {entry_number} contains a starred parameter. "
            "Clean the fit.log before calculating statistics."
        )

    lines = entry.splitlines()
    measurements: dict[ParameterKey, Measurement] = {}
    component_number = 0
    line_index = 0

    while line_index < len(lines):
        model_match = MODEL_PATTERN.match(lines[line_index])

        if model_match is None:
            line_index += 1
            continue

        model_name = model_match.group(1).lower()
        fitted_tokens = extract_number_tokens(model_match.group(2))

        error_line_index = line_index + 1
        while error_line_index < len(lines) and not lines[error_line_index].strip():
            error_line_index += 1

        if error_line_index >= len(lines):
            raise ValueError(
                f"Missing uncertainty line for model '{model_name}' "
                f"in entry {entry_number}."
            )

        fitted_errors = extract_numbers(lines[error_line_index])

        if not fitted_errors:
            raise ValueError(
                f"Could not read uncertainties for model '{model_name}' "
                f"in entry {entry_number}."
            )

        if len(fitted_tokens) < len(fitted_errors):
            raise ValueError(
                f"Model '{model_name}' in entry {entry_number} has "
                f"{len(fitted_tokens)} values but {len(fitted_errors)} errors."
            )

        # The sky line includes two reference coordinates before its three
        # fitted values. Keeping the last N tokens removes those coordinates.
        fitted_tokens = fitted_tokens[-len(fitted_errors) :]
        parameter_names = get_parameter_names(model_name, len(fitted_tokens))

        component_number += 1

        for parameter_index, (parameter_name, token, fitted_error) in enumerate(
            zip(parameter_names, fitted_tokens, fitted_errors),
            start=1,
        ):
            fitted_value, fit_toggle = token
            key = ParameterKey(
                component=component_number,
                model=model_name,
                index=parameter_index,
                name=parameter_name,
            )
            measurements[key] = Measurement(
                value=fitted_value,
                error=abs(fitted_error),
                toggle=fit_toggle,
            )

        line_index = error_line_index + 1

    if not measurements:
        raise ValueError(
            f"No GALFIT model parameters were found in entry {entry_number}."
        )

    return measurements


def circular_mean(values: Sequence[float], period: float = 180.0) -> float:
    """Calculate a circular mean for ellipse position angles."""
    scale = 2.0 * math.pi / period
    sine_sum = sum(math.sin(value * scale) for value in values)
    cosine_sum = sum(math.cos(value * scale) for value in values)
    mean_angle = math.atan2(sine_sum, cosine_sum) / scale

    # Use the representation closest to the first fitted value.
    mean_angle += round((values[0] - mean_angle) / period) * period
    return mean_angle


def circular_deviation(
    value: float,
    mean: float,
    period: float = 180.0,
) -> float:
    """Return the shortest difference between equivalent ellipse angles."""
    half_period = period / 2.0
    return (value - mean + half_period) % period - half_period


def calculate_parameter_statistics(
    parsed_entries: list[dict[ParameterKey, Measurement]],
) -> list[ParameterStatistics]:
    """Calculate mean values and propagated uncertainties."""
    if not parsed_entries:
        return []

    reference_keys = list(parsed_entries[0].keys())
    reference_key_set = set(reference_keys)

    for entry_number, entry in enumerate(parsed_entries[1:], start=2):
        entry_key_set = set(entry.keys())
        if entry_key_set != reference_key_set:
            missing = reference_key_set - entry_key_set
            additional = entry_key_set - reference_key_set
            raise ValueError(
                f"Entry {entry_number} has a different model structure. "
                f"Missing: {sorted(missing, key=str)}; "
                f"additional: {sorted(additional, key=str)}"
            )

    number_of_fits = len(parsed_entries)
    results: list[ParameterStatistics] = []

    for key in reference_keys:
        values = [entry[key].value for entry in parsed_entries]
        errors = [entry[key].error for entry in parsed_entries]
        toggles = [entry[key].toggle for entry in parsed_entries]

        if len(set(toggles)) != 1:
            raise ValueError(
                f"Fit toggle changes between entries for component "
                f"{key.component}, parameter '{key.name}': {toggles}"
            )

        if key.name == "position_angle":
            mean_value = circular_mean(values)
            deviations = [circular_deviation(value, mean_value) for value in values]
            if number_of_fits > 1:
                scatter_standard_deviation = math.sqrt(
                    sum(deviation**2 for deviation in deviations)
                    / (number_of_fits - 1)
                )
            else:
                scatter_standard_deviation = 0.0
        else:
            mean_value = statistics.fmean(values)
            scatter_standard_deviation = (
                statistics.stdev(values) if number_of_fits > 1 else 0.0
            )

        formal_error_on_mean = (
            math.sqrt(sum(error**2 for error in errors)) / number_of_fits
        )
        scatter_error_on_mean = (
            scatter_standard_deviation / math.sqrt(number_of_fits)
            if number_of_fits > 1
            else 0.0
        )
        combined_error = math.hypot(
            formal_error_on_mean,
            scatter_error_on_mean,
        )

        results.append(
            ParameterStatistics(
                key=key,
                toggle=toggles[0],
                number_of_fits=number_of_fits,
                mean=mean_value,
                formal_error_on_mean=formal_error_on_mean,
                scatter_standard_deviation=scatter_standard_deviation,
                scatter_error_on_mean=scatter_error_on_mean,
                combined_error=combined_error,
                minimum=min(values),
                maximum=max(values),
            )
        )

    return results


def write_csv(output_file: Path, results: Sequence[ParameterStatistics]) -> None:
    """Write parameter statistics to CSV."""
    fieldnames = [
        "component",
        "model",
        "parameter",
        "parameter_index",
        "fit_toggle",
        "number_of_fits",
        "mean",
        "combined_error",
        "formal_error_on_mean",
        "scatter_standard_deviation",
        "scatter_error_on_mean",
        "minimum",
        "maximum",
    ]

    with output_file.open("w", encoding="utf-8", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for result in results:
            writer.writerow(
                {
                    "component": result.key.component,
                    "model": result.key.model,
                    "parameter": result.key.name,
                    "parameter_index": result.key.index,
                    "fit_toggle": result.toggle,
                    "number_of_fits": result.number_of_fits,
                    "mean": f"{result.mean:.12g}",
                    "combined_error": f"{result.combined_error:.12g}",
                    "formal_error_on_mean": (f"{result.formal_error_on_mean:.12g}"),
                    "scatter_standard_deviation": (
                        f"{result.scatter_standard_deviation:.12g}"
                    ),
                    "scatter_error_on_mean": (f"{result.scatter_error_on_mean:.12g}"),
                    "minimum": f"{result.minimum:.12g}",
                    "maximum": f"{result.maximum:.12g}",
                }
            )


def format_parameter_value(value: float, parameter_name: str) -> str:
    """Format mean values appropriately for a GALFIT input file."""
    if parameter_name in {"sky", "dsky_dx", "dsky_dy"}:
        return f"{value:.8e}"
    return f"{value:.6f}"


def group_results_by_component(
    results: Sequence[ParameterStatistics],
) -> dict[int, list[ParameterStatistics]]:
    """Group ordered statistics by global GALFIT component number."""
    grouped: dict[int, list[ParameterStatistics]] = {}
    for result in results:
        grouped.setdefault(result.key.component, []).append(result)
    return grouped


def write_component(
    component_number: int,
    component_results: Sequence[ParameterStatistics],
) -> list[str]:
    """Create one GALFIT component block from mean parameters."""
    model = component_results[0].key.model
    line_definitions = MODEL_OUTPUT_LINES.get(model)

    if line_definitions is None:
        raise ValueError(f"No GALFIT output mapping is defined for model '{model}'.")

    by_index = {result.key.index: result for result in component_results}
    lines = [
        f"# Component number: {component_number}",
        f" 0) {model:<22s} # Component type",
    ]

    for line_definition in line_definitions:
        # Standard Sérsic files place the unused parameters 6--8 before
        # axis ratio (9) and position angle (10).
        if line_definition.galfit_number == "9":
            lines.extend(MODEL_EXTRA_LINES.get(model, ()))

        selected = [by_index[index] for index in line_definition.parameter_indices]
        values = " ".join(
            format_parameter_value(result.mean, result.key.name) for result in selected
        )
        toggles = " ".join(str(result.toggle) for result in selected)
        uncertainties = ", ".join(f"{result.combined_error:.6g}" for result in selected)
        lines.append(
            f"{line_definition.galfit_number:>2s}) "
            f"{values} {toggles}  # {line_definition.description}; "
            f"mean error(s): {uncertainties}"
        )

    lines.append(" Z) 0                      # Skip this model? (yes=1, no=0)")
    return lines


def write_galfit_file(
    output_file: Path,
    results: Sequence[ParameterStatistics],
    header: LogHeader,
    source_log: Path,
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
    """Write a complete GALFIT input file using the mean parameters."""
    xmin, xmax, ymin, ymax = header.fit_region
    conv_x, conv_y = convolution_box
    scale_x, scale_y = plate_scale

    lines = [
        f"# Mean-parameter GALFIT file generated from: {source_log}",
        f"# Number of fitted realizations: {results[0].number_of_fits}",
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
        f"K) {scale_x:.6f} {scale_y:<36.6f} # Plate scale (dx dy) [arcsec/pixel]",
        f"O) {display_type:<45s} # Display type",
        f"P) {run_option:<45d} # 0=optimize, 1=model, 2=imgblock, 3=subcomps",
        "",
        "# INITIAL FITTING PARAMETERS",
        "# Values are means from fit.log; uncertainties are included in comments.",
        "# [] in fit.log -> toggle 0; {} -> toggle 2; ordinary value -> toggle 1.",
        "",
    ]

    grouped = group_results_by_component(results)
    for component_number in sorted(grouped):
        lines.extend(write_component(component_number, grouped[component_number]))
        lines.append("")

    lines.append(
        "================================================================================"
    )
    output_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


def print_results(results: Sequence[ParameterStatistics]) -> None:
    """Print a compact statistics summary."""
    current_component: tuple[int, str] | None = None

    for result in results:
        component = (result.key.component, result.key.model)
        if component != current_component:
            current_component = component
            print(f"\nComponent {result.key.component}: {result.key.model}")
            print("-" * 78)

        print(
            f"{result.key.name:24s} = "
            f"{result.mean:14.7g} +/- {result.combined_error:10.4g} "
            f"toggle={result.toggle}"
        )


def validate_headers(headers: Sequence[LogHeader]) -> LogHeader:
    """Verify that all fit.log entries refer to the same image and region."""
    reference = headers[0]
    for entry_number, header in enumerate(headers[1:], start=2):
        if header != reference:
            raise ValueError(
                f"Entry {entry_number} has different input/output image or "
                "fitting-region information."
            )
    return reference


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate mean GALFIT parameters and uncertainties from fit.log, "
            "then write both a CSV table and a GALFIT input file."
        )
    )

    parser.add_argument("fit_log", type=Path, help="GALFIT fit.log file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("fit_statistics.csv"),
        help="Statistics CSV output (default: fit_statistics.csv)",
    )
    parser.add_argument(
        "--galfit-output",
        type=Path,
        default=Path("galfit_mean.init"),
        help="Mean-parameter GALFIT file (default: galfit_mean.init)",
    )

    # Optional overrides for values normally recovered from fit.log.
    parser.add_argument("--input-image", help="Override GALFIT A)")
    parser.add_argument("--output-image", help="Override GALFIT B)")
    parser.add_argument(
        "--fit-region",
        nargs=4,
        type=int,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Override GALFIT H)",
    )

    # Header values absent from fit.log.
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
        default="mask.fits",
        help="GALFIT F) value (default: mask.fits)",
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


def main_meanParamUncer() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not args.fit_log.is_file():
        parser.error(f"Input file does not exist: {args.fit_log}")

    log_text = args.fit_log.read_text(encoding="utf-8")
    entries = split_entries(log_text)

    if not entries:
        parser.error(f"No GALFIT entries were found in {args.fit_log}")

    try:
        parsed_entries = [
            parse_entry(entry, entry_number)
            for entry_number, entry in enumerate(entries, start=1)
        ]
        results = calculate_parameter_statistics(parsed_entries)

        # A, B, and H are read from fit.log unless explicitly overridden.
        if args.fit_region is None:
            headers = [extract_log_header(entry) for entry in entries]
            header = validate_headers(headers)
        else:
            first_header = extract_log_header(entries[0])
            header = LogHeader(
                input_image=first_header.input_image,
                output_image=first_header.output_image,
                fit_region=tuple(args.fit_region),
            )

        header = LogHeader(
            input_image=args.input_image or header.input_image,
            output_image=args.output_image or header.output_image,
            fit_region=header.fit_region,
        )

        write_csv(args.output, results)
        write_galfit_file(
            output_file=args.galfit_output,
            results=results,
            header=header,
            source_log=args.fit_log,
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

    except ValueError as error:
        parser.error(str(error))

    print(f"Entries analyzed : {len(entries)}")
    print_results(results)
    print(f"\nCSV output       : {args.output}")
    print(f"GALFIT output    : {args.galfit_output}")


if __name__ == "__main__":
    main_meanParamUncer()
