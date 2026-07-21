#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


NUMBER_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")

MODEL_PATTERN = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*(.*)$")

STARRED_PARAMETER_PATTERN = re.compile(r"\*[^*\n]+\*")


# Names for common GALFIT model parameters.
# Unknown models are assigned names such as parameter_1, parameter_2, etc.
KNOWN_PARAMETER_NAMES = {
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


@dataclass(frozen=True)
class ParameterKey:
    """Unique identification of a GALFIT model parameter."""

    model: str
    component: int
    index: int
    name: str


@dataclass(frozen=True)
class Measurement:
    """One fitted value and its corresponding GALFIT uncertainty."""

    value: float
    error: float


@dataclass(frozen=True)
class ParameterStatistics:
    """Statistics calculated from all GALFIT entries."""

    key: ParameterKey
    number_of_fits: int
    mean: float
    formal_error_on_mean: float
    scatter_standard_deviation: float
    scatter_error_on_mean: float
    combined_error: float
    minimum: float
    maximum: float


def extract_numbers(text: str) -> list[float]:
    """
    Extract floating-point numbers from a line.

    Characters such as [], {}, (), commas, and asterisks are ignored.
    """
    return [float(number) for number in NUMBER_PATTERN.findall(text)]


def split_entries(log_text: str) -> list[str]:
    """
    Divide fit.log into complete GALFIT entries.

    An entry starts with 'Input image' and continues until the next
    'Input image' line or the end of the file.
    """
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


def get_parameter_names(
    model: str,
    number_of_parameters: int,
) -> tuple[str, ...]:
    """Return descriptive names for the parameters of a GALFIT model."""
    known_names = KNOWN_PARAMETER_NAMES.get(model)

    if known_names is not None and len(known_names) == number_of_parameters:
        return known_names

    return tuple(f"parameter_{index}" for index in range(1, number_of_parameters + 1))


def parse_entry(
    entry: str,
    entry_number: int,
) -> dict[ParameterKey, Measurement]:
    """
    Extract model values and uncertainties from one GALFIT entry.

    Each model line is paired with the uncertainty line immediately
    below it.
    """
    if STARRED_PARAMETER_PATTERN.search(entry):
        raise ValueError(
            f"Entry {entry_number} contains a starred parameter. "
            "Remove problematic entries before calculating statistics."
        )

    lines = entry.splitlines()
    model_occurrences: defaultdict[str, int] = defaultdict(int)

    measurements: dict[ParameterKey, Measurement] = {}

    line_index = 0

    while line_index < len(lines):
        model_match = MODEL_PATTERN.match(lines[line_index])

        if model_match is None:
            line_index += 1
            continue

        model_name = model_match.group(1).lower()
        fitted_values = extract_numbers(model_match.group(2))

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
                f"Could not read uncertainties for model "
                f"'{model_name}' in entry {entry_number}."
            )

        if len(fitted_values) < len(fitted_errors):
            raise ValueError(
                f"Model '{model_name}' in entry {entry_number} "
                f"has {len(fitted_values)} values but "
                f"{len(fitted_errors)} uncertainties."
            )

        # The sky line contains two additional reference coordinates:
        #
        # sky : [x_center, y_center] sky dsky_dx dsky_dy
        #
        # Its uncertainty line contains only the final three parameters.
        # Selecting the last N values removes those reference coordinates.
        fitted_values = fitted_values[-len(fitted_errors) :]

        model_occurrences[model_name] += 1
        component_number = model_occurrences[model_name]

        parameter_names = get_parameter_names(
            model_name,
            len(fitted_values),
        )

        for parameter_index, (parameter_name, fitted_value, fitted_error,) in enumerate(
            zip(
                parameter_names,
                fitted_values,
                fitted_errors,
            ),
            start=1,
        ):
            key = ParameterKey(
                model=model_name,
                component=component_number,
                index=parameter_index,
                name=parameter_name,
            )

            measurements[key] = Measurement(
                value=fitted_value,
                error=abs(fitted_error),
            )

        line_index = error_line_index + 1

    if not measurements:
        raise ValueError(
            f"No GALFIT model parameters were found in entry " f"{entry_number}."
        )

    return measurements


def circular_mean(
    values: list[float],
    period: float = 180.0,
) -> float:
    """
    Calculate a circular mean for quantities with a given period.

    GALFIT position angles describe ellipse orientations, so they have
    a period of 180 degrees rather than 360 degrees.
    """
    scale = 2.0 * math.pi / period

    sine_sum = sum(math.sin(value * scale) for value in values)

    cosine_sum = sum(math.cos(value * scale) for value in values)

    mean_angle = (
        math.atan2(
            sine_sum,
            cosine_sum,
        )
        / scale
    )

    # Select the equivalent representation nearest the first value.
    mean_angle += round((values[0] - mean_angle) / period) * period

    return mean_angle


def circular_deviation(
    value: float,
    mean: float,
    period: float = 180.0,
) -> float:
    """Return the shortest circular difference between two angles."""
    half_period = period / 2.0

    return (value - mean + half_period) % period - half_period


def calculate_parameter_statistics(
    parsed_entries: list[dict[ParameterKey, Measurement]],
) -> list[ParameterStatistics]:
    """Calculate means and propagated uncertainties."""
    if not parsed_entries:
        return []

    reference_keys = list(parsed_entries[0].keys())
    reference_key_set = set(reference_keys)

    for entry_number, entry in enumerate(
        parsed_entries[1:],
        start=2,
    ):
        entry_key_set = set(entry.keys())

        if entry_key_set != reference_key_set:
            missing_parameters = reference_key_set - entry_key_set
            extra_parameters = entry_key_set - reference_key_set

            message = [
                f"Entry {entry_number} does not have the same "
                "model structure as entry 1."
            ]

            if missing_parameters:
                message.append(f"Missing parameters: {missing_parameters}")

            if extra_parameters:
                message.append(f"Additional parameters: {extra_parameters}")

            raise ValueError("\n".join(message))

    number_of_fits = len(parsed_entries)
    statistics_results: list[ParameterStatistics] = []

    for key in reference_keys:
        values = [entry[key].value for entry in parsed_entries]

        errors = [entry[key].error for entry in parsed_entries]

        if key.name == "position_angle":
            mean_value = circular_mean(
                values,
                period=180.0,
            )

            deviations = [
                circular_deviation(
                    value,
                    mean_value,
                    period=180.0,
                )
                for value in values
            ]

            if number_of_fits > 1:
                scatter_standard_deviation = math.sqrt(
                    sum(deviation**2 for deviation in deviations)
                    / (number_of_fits - 1)
                )
            else:
                scatter_standard_deviation = 0.0

        else:
            mean_value = statistics.fmean(values)

            if number_of_fits > 1:
                scatter_standard_deviation = statistics.stdev(values)
            else:
                scatter_standard_deviation = 0.0

        # Propagation of the individual GALFIT uncertainties through
        # an arithmetic mean.
        formal_error_on_mean = (
            math.sqrt(sum(error**2 for error in errors)) / number_of_fits
        )

        # Uncertainty of the mean caused by the variation between
        # different initial-parameter fits.
        if number_of_fits > 1:
            scatter_error_on_mean = scatter_standard_deviation / math.sqrt(
                number_of_fits
            )
        else:
            scatter_error_on_mean = 0.0

        combined_error = math.hypot(
            formal_error_on_mean,
            scatter_error_on_mean,
        )

        statistics_results.append(
            ParameterStatistics(
                key=key,
                number_of_fits=number_of_fits,
                mean=mean_value,
                formal_error_on_mean=formal_error_on_mean,
                scatter_standard_deviation=(scatter_standard_deviation),
                scatter_error_on_mean=scatter_error_on_mean,
                combined_error=combined_error,
                minimum=min(values),
                maximum=max(values),
            )
        )

    return statistics_results


def write_csv(
    output_file: Path,
    results: list[ParameterStatistics],
) -> None:
    """Write the calculated parameter statistics to a CSV file."""
    fieldnames = [
        "model",
        "component",
        "parameter",
        "parameter_index",
        "number_of_fits",
        "mean",
        "combined_error",
        "formal_error_on_mean",
        "scatter_standard_deviation",
        "scatter_error_on_mean",
        "minimum",
        "maximum",
    ]

    with output_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as csv_file:
        writer = csv.DictWriter(
            csv_file,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        for result in results:
            writer.writerow(
                {
                    "model": result.key.model,
                    "component": result.key.component,
                    "parameter": result.key.name,
                    "parameter_index": result.key.index,
                    "number_of_fits": result.number_of_fits,
                    "mean": f"{result.mean:.12g}",
                    "combined_error": (f"{result.combined_error:.12g}"),
                    "formal_error_on_mean": (f"{result.formal_error_on_mean:.12g}"),
                    "scatter_standard_deviation": (
                        f"{result.scatter_standard_deviation:.12g}"
                    ),
                    "scatter_error_on_mean": (f"{result.scatter_error_on_mean:.12g}"),
                    "minimum": f"{result.minimum:.12g}",
                    "maximum": f"{result.maximum:.12g}",
                }
            )


def print_results(
    results: list[ParameterStatistics],
) -> None:
    """Print a readable summary to the terminal."""
    current_component: tuple[str, int] | None = None

    for result in results:
        component = (
            result.key.model,
            result.key.component,
        )

        if component != current_component:
            current_component = component

            print(f"\n{result.key.model} " f"component {result.key.component}")

            print("-" * 72)

        print(
            f"{result.key.name:24s} = "
            f"{result.mean:14.7g} +/- "
            f"{result.combined_error:10.4g}"
            f"   formal={result.formal_error_on_mean:.4g}"
            f"   scatter={result.scatter_standard_deviation:.4g}"
        )


def main_meanParamUncer() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate mean GALFIT parameters and propagated "
            "uncertainties from multiple fit.log entries."
        )
    )

    parser.add_argument(
        "fit_log",
        type=Path,
        help="GALFIT fit.log file",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("fit_statistics.csv"),
        help=("Output CSV file " "(default: fit_statistics.csv)"),
    )

    args = parser.parse_args()

    if not args.fit_log.is_file():
        parser.error(f"Input file does not exist: {args.fit_log}")

    log_text = args.fit_log.read_text(
        encoding="utf-8",
    )

    entries = split_entries(log_text)

    if not entries:
        parser.error(f"No GALFIT entries were found in {args.fit_log}")

    try:
        parsed_entries = [
            parse_entry(
                entry,
                entry_number,
            )
            for entry_number, entry in enumerate(
                entries,
                start=1,
            )
        ]

        results = calculate_parameter_statistics(parsed_entries)

    except ValueError as error:
        parser.error(str(error))

    write_csv(
        args.output,
        results,
    )

    print(f"Entries analyzed: {len(entries)}")
    print_results(results)
    print(f"\nCSV output: {args.output}")


if __name__ == "__main__":
    main_meanParamUncer()
