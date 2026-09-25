#!/usr/bin/env python3
"""
Prepare TOGA gene gain/loss data for trait correlation using BayesTraits.

This script merges a gene-state table with a trait table, creates a BayesTraits
input file, and optionally runs BayesTraits if provided.

Example:
    python toga_trait_analysis.py \
        --gene-states gene_states.tsv \
        --traits traits.tsv \
        --trait-column max_longevity_y \
        --tree /path/to/tree.newick \
        --output-dir results/
"""

import argparse
import csv
import statistics
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple


def sniff_delimiter(file_path: Path) -> str:
    with open(file_path, "r", newline="") as f:
        sample = f.read(2048)
    sniffer = csv.Sniffer()
    try:
        dialect = sniffer.sniff(sample, delimiters=[",", "\t", ";"])
        return dialect.delimiter
    except csv.Error:
        return "\t"


def load_table(file_path: Path) -> Tuple[List[dict], str]:
    delimiter = sniff_delimiter(file_path)
    with open(file_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        rows = [row for row in reader]
    return rows, delimiter


def normalize_species(name: str) -> str:
    return name.strip().replace(" ", "_")


def parse_state(value: str) -> int:
    if value is None:
        return 0
    val = str(value).strip().lower()
    if val in {"1", "gain", "gained", "dup", "duplication", "present", "yes"}:
        return 1
    if val in {"0", "loss", "lost", "absent", "no"}:
        return 0
    try:
        return 1 if float(val) > 0 else 0
    except ValueError:
        return 0


def build_state_map(rows: List[dict], species_column: str, state_column: str) -> Dict[str, int]:
    state_map = {}
    for row in rows:
        species = normalize_species(row.get(species_column, ""))
        if not species:
            continue
        state_map[species] = parse_state(row.get(state_column, ""))
    return state_map


def build_trait_map(rows: List[dict], species_column: str, trait_column: str) -> Dict[str, float]:
    trait_map = {}
    for row in rows:
        species = normalize_species(row.get(species_column, ""))
        if not species:
            continue
        value = row.get(trait_column, "")
        try:
            trait_map[species] = float(value)
        except (TypeError, ValueError):
            continue
    return trait_map


def write_csv(output_file: Path, rows: List[dict]):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_bayestraits_input(output_file: Path, rows: List[dict], columns: List[str]):
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in rows:
            writer.writerow([row[col] for col in columns])


def run_bayestraits(bayestraits_path: Path, tree_file: Path, data_file: Path,
                    command_file: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    stdout_file = output_dir / "bayestraits_stdout.txt"
    stderr_file = output_dir / "bayestraits_stderr.txt"

    with open(command_file, "r") as cmd_in, \
            open(stdout_file, "w") as out, \
            open(stderr_file, "w") as err:
        subprocess.run(
            [str(bayestraits_path), str(tree_file), str(data_file)],
            stdin=cmd_in,
            stdout=out,
            stderr=err,
            check=False,
            text=True,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Prepare TOGA gene-state data for BayesTraits analysis"
    )
    parser.add_argument("--gene-states", required=True, type=Path,
                        help="CSV/TSV with species + gene state (gain/loss)")
    parser.add_argument("--traits", required=True, type=Path,
                        help="CSV/TSV with species traits (e.g. max_longevity_y)")
    parser.add_argument("--trait-column", default="max_longevity_y",
                        help="Trait column name in traits file")
    parser.add_argument("--species-column", default="species",
                        help="Species column name (both files)")
    parser.add_argument("--state-column", default="state",
                        help="State column name in gene-states file")
    parser.add_argument("--discretize-trait", action="store_true",
                        help="Discretize trait into 0/1 using a threshold")
    parser.add_argument("--trait-threshold", type=float, default=None,
                        help="Threshold for discretizing trait (default: median)")
    parser.add_argument("--tree", type=Path, default=None,
                        help="Tree file for BayesTraits (Newick)")
    parser.add_argument("--bayestraits-path", type=Path, default=None,
                        help="Path to BayesTraits executable (optional)")
    parser.add_argument("--bayestraits-commands", type=Path, default=None,
                        help="Command file for BayesTraits batch mode")
    parser.add_argument("--output-dir", type=Path, default=Path("output/toga_trait_analysis"),
                        help="Output directory (default: output/toga_trait_analysis)")

    args = parser.parse_args()

    gene_rows, _ = load_table(args.gene_states)
    trait_rows, _ = load_table(args.traits)

    state_map = build_state_map(gene_rows, args.species_column, args.state_column)
    trait_map = build_trait_map(trait_rows, args.species_column, args.trait_column)

    shared_species = sorted(set(state_map.keys()) & set(trait_map.keys()))
    if not shared_species:
        raise SystemExit("No overlapping species between gene states and traits.")

    traits = [trait_map[sp] for sp in shared_species]
    threshold = args.trait_threshold
    if args.discretize_trait and threshold is None:
        threshold = statistics.median(traits)

    merged_rows = []
    for sp in shared_species:
        trait_value = trait_map[sp]
        trait_binary = 1 if args.discretize_trait and trait_value >= threshold else 0
        merged_rows.append({
            "species": sp,
            "gene_state": state_map[sp],
            "trait_value": trait_value,
            "trait_binary": trait_binary,
        })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    merged_file = args.output_dir / "merged_trait_data.csv"
    write_csv(merged_file, merged_rows)

    if args.discretize_trait:
        bayestraits_columns = ["species", "gene_state", "trait_binary"]
    else:
        bayestraits_columns = ["species", "gene_state", "trait_value"]

    bayestraits_file = args.output_dir / "bayestraits_input.txt"
    write_bayestraits_input(bayestraits_file, merged_rows, bayestraits_columns)

    example_cmd_file = args.output_dir / "bayestraits_commands.example.txt"
    if not example_cmd_file.exists():
        example_cmd_file.write_text(
            "# Example BayesTraits batch commands\n"
            "# Replace with the model you want (Discrete/Continuous).\n"
            "# For Discrete binary traits:\n"
            "# 1\n"
            "# Independent\n"
            "# Run\n"
        )

    if args.bayestraits_path and args.tree and args.bayestraits_commands:
        run_bayestraits(args.bayestraits_path, args.tree, bayestraits_file,
                        args.bayestraits_commands, args.output_dir)


if __name__ == "__main__":
    main()
