#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


SEPARATOR_PATTERN = re.compile(r"^-{20,}\s*$")
ASTERISK_PARAMETER_PATTERN = re.compile(r"\*[^*\n]+\*")


def extract_entries(log_text: str) -> tuple[str, list[str]]:
    """Extract individual GALFIT entries from a fit.log file."""
    lines = log_text.splitlines()

    separator = "-" * 77
    entries = []
    current_entry = []
    reading_entry = False

    for line in lines:
        if SEPARATOR_PATTERN.match(line):
            separator = line.rstrip()

            if reading_entry and current_entry:
                entries.append("\n".join(current_entry).strip())
                current_entry = []
                reading_entry = False

            continue

        if line.lstrip().startswith("Input image"):
            if current_entry:
                entries.append("\n".join(current_entry).strip())

            current_entry = [line]
            reading_entry = True

        elif reading_entry:
            current_entry.append(line)

    if current_entry:
        entries.append("\n".join(current_entry).strip())

    return separator, entries


def remove_flagged_entries(entries: list[str]) -> tuple[list[str], list[str]]:
    """Separate valid entries from entries containing starred parameters."""
    valid_entries = []
    removed_entries = []

    for entry in entries:
        if ASTERISK_PARAMETER_PATTERN.search(entry):
            removed_entries.append(entry)
        else:
            valid_entries.append(entry)

    return valid_entries, removed_entries


def format_entries(entries: list[str], separator: str) -> str:
    """Reconstruct the GALFIT log using the original separator."""
    if not entries:
        return ""

    between_entries = f"\n\n{separator}\n{separator}\n\n"
    output = f"{separator}\n\n"
    output += between_entries.join(entries)
    output += f"\n\n{separator}\n"

    return output


def clean_fit_log(input_file: Path, output_file: Path) -> None:
    """Remove complete GALFIT entries containing starred parameters."""
    log_text = input_file.read_text(encoding="utf-8")

    separator, entries = extract_entries(log_text)
    valid_entries, removed_entries = remove_flagged_entries(entries)

    cleaned_text = format_entries(valid_entries, separator)
    output_file.write_text(cleaned_text, encoding="utf-8")

    print(f"Total entries   : {len(entries)}")
    print(f"Entries retained: {len(valid_entries)}")
    print(f"Entries removed : {len(removed_entries)}")
    print(f"Output file     : {output_file}")


def mainCleanFitlog() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Remove complete GALFIT fit.log entries containing parameters "
            "enclosed by asterisks."
        )
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="Input GALFIT fit.log file",
    )

    parser.add_argument(
        "output_file",
        type=Path,
        help="Output cleaned log file",
    )

    args = parser.parse_args()

    if not args.input_file.is_file():
        parser.error(f"Input file does not exist: {args.input_file}")

    clean_fit_log(args.input_file, args.output_file)


if __name__ == "__main__":
    mainCleanFitlog()
