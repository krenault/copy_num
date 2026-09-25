#!/usr/bin/env python3
"""
Batch validate multiple genes using NCBI BLAST RBH.

Reads a TSV file with columns: gene_symbol, target_species
and runs validation for each. Outputs CSV files plus a CSV summary.

Usage:
    python batch_validate.py --input genes_to_validate.tsv --output-dir results/
"""

import argparse
import csv
import time
from pathlib import Path

from validate_gene import validate_gene_single, OUTPUT_DIR


def write_result_csv(result: dict, output_file: Path) -> None:
    fieldnames = [
        "gene_symbol",
        "target_species",
        "status",
        "reference_source",
        "reference_id",
        "reference_length",
        "n_hits",
        "ecnc",
        "duplicated_call",
        "top_hit_accession",
        "top_hit_identity",
        "top_hit_evalue",
        "top_hit_definition",
        "reciprocal_status",
        "reciprocal_top_hit_accession",
        "reciprocal_top_hit_identity",
        "reciprocal_top_hit_evalue",
        "reciprocal_top_hit_definition",
    ]
    row = {key: result.get(key, "") for key in fieldnames}
    if "reciprocal_details" in result:
        row["reciprocal_top_hit_accession"] = result["reciprocal_details"].get("top_hit_accession", "")
        row["reciprocal_top_hit_identity"] = result["reciprocal_details"].get("top_hit_identity", "")
        row["reciprocal_top_hit_evalue"] = result["reciprocal_details"].get("top_hit_evalue", "")
        row["reciprocal_top_hit_definition"] = result["reciprocal_details"].get("top_hit_definition", "")

    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)


def load_gene_list(input_file: Path) -> list:
    """Load genes to validate from TSV file."""
    genes = []

    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genes.append({
                "gene_symbol": row.get("gene_symbol", row.get("gene", "")),
                "target_species": row.get("target_species", row.get("species", "")),
            })

    return genes


def batch_validate(gene_list: list, output_dir: Path, delay: int = 60) -> list:
    """Validate multiple genes with delay between each."""

    all_results = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for i, gene_info in enumerate(gene_list):
        print(f"\n{'#' * 70}")
        print(f"# Processing {i+1}/{len(gene_list)}: {gene_info['gene_symbol']}")
        print(f"{'#' * 70}")

        results = validate_gene_single(
            gene_symbol=gene_info["gene_symbol"],
            target_species=gene_info["target_species"],
            reference_source="auto",
            use_protein=False,
        )

        # Save individual results as CSV
        target_clean = gene_info["target_species"].replace(" ", "_").lower()
        output_file = output_dir / f"{gene_info['gene_symbol']}_{target_clean}_validation.csv"
        write_result_csv(results, output_file)

        all_results.append({
            "gene_symbol": gene_info["gene_symbol"],
            "target_species": gene_info["target_species"],
            "status": results.get("status", ""),
            "top_hit_accession": results.get("top_hit_accession", ""),
            "reciprocal_status": results.get("reciprocal_status", ""),
        })

        # Wait between genes to avoid rate limiting
        if i < len(gene_list) - 1:
            print(f"\nWaiting {delay}s before next gene...")
            time.sleep(delay)

    # Save summary CSV
    summary_file = output_dir / "batch_validation_summary.csv"
    if all_results:
        fieldnames = list(all_results[0].keys())
        with open(summary_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)

    print(f"\n{'=' * 70}")
    print("BATCH VALIDATION COMPLETE")
    print(f"{'=' * 70}")
    print(f"Processed {len(gene_list)} genes")
    print(f"Results saved to: {output_dir}")

    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Batch validate multiple genes using NCBI BLAST"
    )

    parser.add_argument("--input", "-i", type=Path, required=True,
                        help="TSV file with genes to validate (columns: gene_symbol, target_species)")
    parser.add_argument("--output-dir", "-o", type=Path, default=OUTPUT_DIR,
                        help="Output directory for results")
    parser.add_argument("--delay", "-d", type=int, default=60,
                        help="Delay in seconds between genes (default: 60)")

    args = parser.parse_args()

    # Load gene list
    gene_list = load_gene_list(args.input)
    print(f"Loaded {len(gene_list)} genes to validate")

    # Run batch validation
    batch_validate(gene_list, args.output_dir, args.delay)


if __name__ == "__main__":
    main()
