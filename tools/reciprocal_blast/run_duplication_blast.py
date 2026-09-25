#!/usr/bin/env python3
"""
Run BLAST validation on top genes from duplication_interest_scores.tsv.

Usage:
    python run_duplication_blast.py --n-genes 20 --focal-species "Myotis davidii"
    python run_duplication_blast.py --n-genes 10 --focal-species "Sciurus carolinensis"
"""

import argparse
import csv
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
INPUT_FILE = PROJECT_ROOT / "output" / "duplication_dollo" / "duplication_interest_scores.tsv"
OUTPUT_DIR = PROJECT_ROOT / "output" / "copy_num_BLAT_check"


def load_top_genes(input_file: Path, n_genes: int, order_filter: str = None) -> list:
    """Load top N genes from interest scores file."""
    genes = []
    with open(input_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if order_filter and row.get("primary_order", "") != order_filter:
                continue
            genes.append({
                "gene_symbol": row["gene_symbol"],
                "gene_id": row["gene_id"],
                "interest_score": float(row.get("interest_score", 0)),
                "primary_order": row.get("primary_order", ""),
                "n_species_dup": int(row.get("n_species_dup", 0)),
            })
            if len(genes) >= n_genes:
                break
    return genes


def main():
    parser = argparse.ArgumentParser(description="Run BLAST validation on top duplication genes")
    parser.add_argument("--n-genes", type=int, default=20, help="Number of top genes to validate")
    parser.add_argument("--focal-species", required=True, help="Focal species for neighbor search")
    parser.add_argument("--n-neighbors", type=int, default=5, help="Number of neighbors to check")
    parser.add_argument("--order-filter", default=None, help="Filter by primary_order (e.g., Chiroptera)")
    parser.add_argument("--input-file", type=Path, default=INPUT_FILE, help="Input TSV file")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    args = parser.parse_args()

    genes = load_top_genes(args.input_file, args.n_genes, args.order_filter)
    print(f"Loaded {len(genes)} genes to validate")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for i, gene in enumerate(genes):
        print(f"\n[{i+1}/{len(genes)}] {gene['gene_symbol']} (score={gene['interest_score']}, order={gene['primary_order']})")

        cmd = [
            "python3", str(SCRIPT_DIR / "validate_gene.py"),
            "--gene", gene["gene_symbol"],
            "--focal-species", args.focal_species,
            "--n-neighbors", str(args.n_neighbors),
            "--require-ncbi",
            "--protein",
        ]

        if args.dry_run:
            print("  Would run:", " ".join(cmd))
        else:
            import subprocess
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                if result.returncode == 0:
                    print(f"  Done: {gene['gene_symbol']}")
                else:
                    print(f"  Error: {result.stderr[:200]}")
            except subprocess.TimeoutExpired:
                print(f"  Timeout for {gene['gene_symbol']}")
            except Exception as e:
                print(f"  Exception: {e}")

            # Rate limiting
            time.sleep(2)

    print(f"\nResults saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
