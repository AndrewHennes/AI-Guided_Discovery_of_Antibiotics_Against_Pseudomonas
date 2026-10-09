#!/usr/bin/env python3
r"""Append efflux hit probabilities to a CSV containing SMILES.

Requires Python 3.9+ and RDKit. No import of the original pipeline is needed.
The default input column is "smiles"; use --smiles-column to change it.
The output preserves input columns and adds "efflux_rules_score" = P(hit).
Missing/invalid SMILES and unavailable descriptors produce blank scores.

Example (replace the input path with your CSV's absolute path):
    python /Users/asselism/Documents/Codex/2026-09-22/ple/outputs/efflux_rules_score.py \
        /absolute/path/compounds.csv --smiles-column SMILES

Without --output, results are saved beside this script as
<input_stem>_efflux_scores.csv, regardless of the current working directory.
An existing output file is replaced; the input file cannot be the output.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import sys
import tempfile
from typing import Optional, Sequence


OUTPUT_DIRECTORY = Path(__file__).resolve().parent
SCORE_COLUMN = "efflux_rules_score"


def posterior_probability_from_multiple_gaussians(
    sample: Sequence[float],
    standard_deviations: Sequence[Sequence[float]],
    means: Sequence[Sequence[float]],
    priors: Sequence[Sequence[float]],
) -> list[float]:
    """Return normalized class probabilities for independent Gaussian dimensions.

    means, standard_deviations, and priors have shape (dimensions, classes).
    Only priors[0] is used, so the class prior is applied once per sample.
    Calculations use log probabilities to avoid likelihood underflow.
    Invalid dimensions or nonfinite parameters raise ValueError.
    """
    number_of_dimensions = len(sample)
    if number_of_dimensions == 0 or len(priors) == 0 or len(priors[0]) == 0:
        raise ValueError("At least one dimension and one class are required.")
    number_of_classes = len(priors[0])
    for name, matrix in (
        ("means", means),
        ("standard_deviations", standard_deviations),
        ("priors", priors),
    ):
        if len(matrix) != number_of_dimensions or any(
            len(row) != number_of_classes for row in matrix
        ):
            raise ValueError(f"{name} must have shape (dimensions, classes).")
    if not all(math.isfinite(value) for value in sample):
        raise ValueError("Sample values must be finite.")
    if not all(math.isfinite(value) for row in means for value in row):
        raise ValueError("Means must be finite.")
    if not all(
        math.isfinite(value) and value > 0
        for row in standard_deviations
        for value in row
    ):
        raise ValueError("Standard deviations must be finite and positive.")
    if not all(math.isfinite(p) and p >= 0 for p in priors[0]) or not any(
        p > 0 for p in priors[0]
    ):
        raise ValueError("Priors must be finite, nonnegative, and not all zero.")

    base_log_probability = -0.5 * number_of_dimensions * math.log(2 * math.pi)
    log_probabilities = [
        math.log(prior) + base_log_probability if prior > 0 else -math.inf
        for prior in priors[0]
    ]
    for value, dimension_means, dimension_standard_deviations in zip(
        sample, means, standard_deviations
    ):
        for class_index in range(number_of_classes):
            standard_deviation = dimension_standard_deviations[class_index]
            standardized_error = (value - dimension_means[class_index]) / standard_deviation
            log_probabilities[class_index] += (
                -math.log(standard_deviation) - 0.5 * standardized_error * standardized_error
            )
    maximum_log_probability = max(log_probabilities)
    if not math.isfinite(maximum_log_probability):
        raise ValueError("Sample is outside the supported floating-point range.")
    unnormalized_probabilities = [
        math.exp(value - maximum_log_probability) for value in log_probabilities
    ]
    normalization = sum(unnormalized_probabilities)
    return [probability / normalization for probability in unnormalized_probabilities]


def apply_efflux_rules_score(
    ipc: float, logp: float, number_aromatics: float
) -> list[float]:
    """Return [P(hit), P(nonhit)]; ipc must already be log10(Ipc)."""
    means = [[4.14, 5.89], [2.43, 3.75], [1.63, 2.55]]
    standard_deviations = [[2.07, 1.36], [1.22, 1.4], [1.4, 0.98]]
    prior_hit = 2.5 / 10
    prior_nonhit = 1 - prior_hit
    priors = [[prior_hit, prior_nonhit] for _ in range(3)]
    return posterior_probability_from_multiple_gaussians(
        [ipc, logp, number_aromatics], standard_deviations, means, priors
    )


def apply_efflux_rules_score_to_smiles(smiles: Optional[str]) -> Optional[float]:
    """Return P(hit) or None when a molecule cannot be scored.

    Match the pipeline descriptors: log10(Descriptors.Ipc(mol)),
    Descriptors.MolLogP(mol), and CalcNumAromaticRings(mol).
    Ipc uses RDKit's default avg=False. No molecular standardization is applied.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None

    # Import lazily so numeric scoring and --help also work without RDKit.
    from rdkit import Chem, rdBase
    from rdkit.Chem import Descriptors, rdMolDescriptors

    try:
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(smiles.strip())
            if molecule is None or molecule.GetNumAtoms() == 0:
                return None
            ipc_value = float(Descriptors.Ipc(molecule))
            if not math.isfinite(ipc_value) or ipc_value <= 0:
                return None
            sample = [
                math.log10(ipc_value),
                float(Descriptors.MolLogP(molecule)),
                int(rdMolDescriptors.CalcNumAromaticRings(molecule)),
            ]
        if not all(math.isfinite(value) for value in sample):
            return None
        return apply_efflux_rules_score(*sample)[0]
    except (ValueError, RuntimeError, ArithmeticError):
        return None


def score_csv(
    input_path: Path,
    output_path: Optional[Path] = None,
    smiles_column: str = "smiles",
) -> tuple[Path, int, int]:
    """Write a scored CSV and return (absolute output path, rows, unscored rows).

    Stream rows without converting other columns' values or dropping rows.
    Commit the output only after the complete input has been processed.
    """
    input_path = Path(input_path).expanduser().resolve()
    output_path = (
        OUTPUT_DIRECTORY / f"{input_path.stem}_efflux_scores.csv"
        if output_path is None
        else Path(output_path).expanduser().resolve()
    )
    if input_path == output_path or (
        output_path.exists() and input_path.samefile(output_path)
    ):
        raise ValueError("Input and output must be different files.")

    # Fail before writing if the required dependency cannot be loaded.
    try:
        from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("RDKit is required. Run with a Python environment containing RDKit.") from exc

    with input_path.open("r", encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source, strict=True)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError("The input CSV has no header.")
        if len(set(fieldnames)) != len(fieldnames):
            raise ValueError("The input CSV contains duplicate column names.")
        if smiles_column not in fieldnames:
            raise ValueError(
                f"SMILES column {smiles_column!r} was not found. "
                f"Available columns: {', '.join(fieldnames)}. Use --smiles-column."
            )
        if SCORE_COLUMN in fieldnames:
            raise ValueError(f"The input already contains the output column {SCORE_COLUMN!r}.")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = None
        total = unscored = 0
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", newline="", delete=False,
                dir=str(output_path.parent), prefix=f".{output_path.name}.", suffix=".tmp",
            ) as destination:
                temporary_path = Path(destination.name).resolve()
                writer = csv.DictWriter(destination, fieldnames=[*fieldnames, SCORE_COLUMN])
                writer.writeheader()
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError(
                            f"Malformed CSV record ending at line {reader.line_num}: "
                            "the number of fields differs from the header."
                        )
                    score = apply_efflux_rules_score_to_smiles(row[smiles_column])
                    row[SCORE_COLUMN] = "" if score is None else score
                    writer.writerow(row)
                    total += 1
                    unscored += score is None
            os.replace(temporary_path, output_path)
        finally:
            if temporary_path is not None:
                temporary_path.unlink(missing_ok=True)
    return output_path, total, unscored


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_csv", type=Path, help="Path to a UTF-8 CSV with a header row.")
    parser.add_argument(
        "--smiles-column", default="smiles", metavar="NAME",
        help="SMILES column name (case-sensitive; default: smiles).",
    )
    parser.add_argument(
        "-o", "--output", type=Path, metavar="CSV",
        help="Output CSV path (default: <script directory>/<input stem>_efflux_scores.csv).",
    )
    args = parser.parse_args()
    try:
        output_path, total, unscored = score_csv(args.input_csv, args.output, args.smiles_column)
    except (OSError, ValueError, RuntimeError, csv.Error) as exc:
        parser.exit(2, f"Error: {exc}\n")
    print(f"Saved {total} rows to {output_path}")
    print(f"Scored: {total - unscored}; blank scores: {unscored}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
