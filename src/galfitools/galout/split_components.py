#!/usr/bin/env python3
"""Split a GALFIT file into one file per non-sky component."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


COMPONENT_RE = re.compile(
    r"^\s*#\s*Component\s+number\s*:\s*(\d+)\s*$",
    re.IGNORECASE,
)

TYPE_RE = re.compile(
    r"^\s*0\)\s*(\S+)",
    re.IGNORECASE,
)

FINAL_SEPARATOR_RE = re.compile(r"^\s*=+\s*$")


def split_galfit_file(
    galfit_file: Path,
    output_dir: Path | None = None,
) -> list[Path]:
    """Split a GALFIT file into several files, one per non-sky component.

    Parameters
    ----------
    galfit_file : pathlib.Path
        Input GALFIT file.
    output_dir : pathlib.Path or None, optional
        Directory where the output files are written. If `None`, the files are
        written in the same directory as `galfit_file`.

    Returns
    -------
    list of pathlib.Path
        Paths of the generated files.

    Notes
    -----
    The script keeps the original GALFIT header, defined as all lines before
    the first ``# Component number:`` block. The sky component is ignored.
    """
    galfit_file = galfit_file.expanduser().resolve()

    if output_dir is None:
        output_dir = galfit_file.parent
    else:
        output_dir = output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

    lines = galfit_file.read_text().splitlines(keepends=True)

    component_starts: list[int] = []
    component_numbers: list[int] = []

    for i, line in enumerate(lines):
        match = COMPONENT_RE.match(line)
        if match:
            component_starts.append(i)
            component_numbers.append(int(match.group(1)))

    if not component_starts:
        raise ValueError(f"No GALFIT component blocks were found in {galfit_file}")

    header = lines[: component_starts[0]]
    generated_files: list[Path] = []

    base_name = galfit_file.stem
    extension = galfit_file.suffix

    for idx, start in enumerate(component_starts):
        end = (
            component_starts[idx + 1] if idx + 1 < len(component_starts) else len(lines)
        )

        block = lines[start:end]

        # Remove empty lines and the final GALFIT separator from the last block.
        while block and (
            block[-1].strip() == "" or FINAL_SEPARATOR_RE.match(block[-1])
        ):
            block = block[:-1]

        component_type = get_component_type(block)

        if component_type.lower() == "sky":
            continue

        component_number = component_numbers[idx]
        output_name = f"{base_name}-{component_type}-nc{component_number}{extension}"
        output_path = output_dir / output_name

        output_lines = []
        output_lines.extend(header)
        output_lines.extend(block)

        if output_lines and not output_lines[-1].endswith("\n"):
            output_lines[-1] += "\n"

        output_lines.append("\n")
        output_lines.append("=" * 80 + "\n")

        output_path.write_text("".join(output_lines))
        generated_files.append(output_path)

    return generated_files


def get_component_type(component_block: list[str]) -> str:
    """Return the GALFIT function name for a component block."""
    for line in component_block:
        match = TYPE_RE.match(line)
        if match:
            return match.group(1).strip()

    raise ValueError("A component block does not contain a valid '0)' type line")


def mainSplitComp() -> None:
    parser = argparse.ArgumentParser(
        description="Split a GALFIT file into one file per non-sky component."
    )

    parser.add_argument(
        "galfit_file",
        help="Input GALFIT file, for example galfit.02",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory. Default: same directory as the input file.",
    )

    args = parser.parse_args()

    output_files = split_galfit_file(
        Path(args.galfit_file),
        Path(args.output_dir) if args.output_dir is not None else None,
    )

    for output_file in output_files:
        print(output_file)


if __name__ == "__main__":
    mainSplitComp()
