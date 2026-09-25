#!/usr/bin/env python3
"""
Fetch genome assembly quality metrics from NCBI for vertebrates.
Creates a CSV of high-quality genomes suitable for copy number analysis.
"""

import requests
import csv
import time
import argparse
from pathlib import Path
from typing import Optional
import xml.etree.ElementTree as ET

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_DELAY = 0.35  # Rate limit: ~3 requests/sec


def fetch_assembly_summaries(
    taxon: str = "Vertebrata",
    min_n50: int = 50_000_000,
    refseq_only: bool = True,
    limit: Optional[int] = None,
    api_key: Optional[str] = None
) -> list[dict]:
    """
    Fetch assembly summaries from NCBI Datasets API.

    Args:
        taxon: Taxonomic group to search
        min_n50: Minimum scaffold N50 in bp
        refseq_only: Only include RefSeq assemblies
        limit: Max assemblies to fetch (None = all)
        api_key: NCBI API key for higher rate limits

    Returns:
        List of assembly dictionaries with quality metrics
    """
    session = requests.Session()
    if api_key:
        session.headers["api-key"] = api_key

    assemblies = []
    page_token = None
    page_num = 0
    total_scanned = 0

    print(f"Fetching {taxon} assemblies with N50 >= {min_n50/1e6:.0f} Mb...")
    print(f"  RefSeq only: {refseq_only}")

    while True:
        page_num += 1

        # Build request URL - simpler params, filter client-side
        url = f"{DATASETS_API}/genome/taxon/{taxon}/dataset_report"
        params = {"page_size": 1000}

        if refseq_only:
            params["filters.reference_only"] = "false"  # Get all RefSeq, not just reference
            params["filters.exclude_paired_reports"] = "true"

        if page_token:
            params["page_token"] = page_token

        try:
            time.sleep(NCBI_DELAY)
            resp = session.get(url, params=params, timeout=120)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"  Error on page {page_num}: {e}")
            break

        reports = data.get("reports", [])
        if not reports:
            break

        for report in reports:
            total_scanned += 1

            # Filter by source
            source = report.get("source_database", "")
            if refseq_only and "REFSEQ" not in source:
                continue

            assembly_info = report.get("assembly_info", {})
            assembly_stats = report.get("assembly_stats", {})
            organism = report.get("organism", {})

            # Filter by assembly level
            level = assembly_info.get("assembly_level", "")
            if level not in ("Chromosome", "Complete Genome"):
                continue

            # Get scaffold N50
            scaffold_n50 = assembly_stats.get("scaffold_n50", 0)
            if scaffold_n50 < min_n50:
                continue

            # Extract key metrics
            entry = {
                "accession": report.get("accession", ""),
                "organism_name": organism.get("organism_name", ""),
                "species_taxid": organism.get("tax_id", ""),
                "common_name": organism.get("common_name", ""),

                # Taxonomy
                "infraclass": "",
                "order": "",
                "family": "",
                "genus": "",

                # Assembly info
                "assembly_name": assembly_info.get("assembly_name", ""),
                "assembly_level": level,
                "assembly_type": assembly_info.get("assembly_type", ""),
                "submission_date": assembly_info.get("release_date", ""),
                "refseq_category": assembly_info.get("refseq_category", ""),

                # Quality metrics
                "scaffold_n50": scaffold_n50,
                "contig_n50": assembly_stats.get("contig_n50", 0),
                "total_sequence_length": assembly_stats.get("total_sequence_length", 0),
                "scaffold_count": assembly_stats.get("number_of_scaffolds", 0),
                "contig_count": assembly_stats.get("number_of_contigs", 0),
                "gc_percent": assembly_stats.get("gc_percent", 0),
                "genome_coverage": assembly_stats.get("genome_coverage", ""),

                # Annotation
                "annotation_name": report.get("annotation_info", {}).get("name", ""),
                "gene_count": report.get("annotation_info", {}).get("stats", {}).get("gene_counts", {}).get("total", 0),
                "protein_coding_count": report.get("annotation_info", {}).get("stats", {}).get("gene_counts", {}).get("protein_coding", 0),
            }

            assemblies.append(entry)

            if limit and len(assemblies) >= limit:
                break

        print(f"  Page {page_num}: scanned {total_scanned}, kept {len(assemblies)} high-quality...", flush=True)

        if limit and len(assemblies) >= limit:
            break

        # Check for next page
        page_token = data.get("next_page_token")
        if not page_token:
            break

    print(f"  Done: {len(assemblies)} assemblies from {total_scanned} total")
    return assemblies


def fetch_taxonomy_info(taxids: list[int], session: requests.Session) -> dict:
    """Fetch taxonomy lineage for a batch of taxids."""
    taxonomy_map = {}

    # Process in batches of 200
    batch_size = 200
    for i in range(0, len(taxids), batch_size):
        batch = taxids[i:i + batch_size]
        ids_str = ",".join(str(t) for t in batch)

        url = f"{NCBI_BASE}/efetch.fcgi"
        params = {
            "db": "taxonomy",
            "id": ids_str,
            "retmode": "xml"
        }

        try:
            time.sleep(NCBI_DELAY)
            resp = session.get(url, params=params, timeout=60)
            root = ET.fromstring(resp.content)

            for taxon in root.findall(".//Taxon"):
                taxid = taxon.findtext("TaxId")
                lineage = {}

                for lin_taxon in taxon.findall(".//LineageEx/Taxon"):
                    rank = lin_taxon.findtext("Rank", "")
                    name = lin_taxon.findtext("ScientificName", "")
                    if rank in ("infraclass", "order", "family", "genus", "superorder", "class"):
                        lineage[rank] = name

                taxonomy_map[int(taxid)] = lineage

        except Exception as e:
            print(f"  Taxonomy fetch error: {e}")
            continue

    return taxonomy_map


def enrich_with_taxonomy(assemblies: list[dict]) -> list[dict]:
    """Add taxonomic lineage information to assemblies."""
    print("Fetching taxonomy information...")

    session = requests.Session()
    taxids = list(set(a["species_taxid"] for a in assemblies if a["species_taxid"]))

    taxonomy_map = fetch_taxonomy_info(taxids, session)

    for assembly in assemblies:
        taxid = assembly["species_taxid"]
        if taxid in taxonomy_map:
            lineage = taxonomy_map[taxid]
            assembly["infraclass"] = lineage.get("infraclass", lineage.get("class", ""))
            assembly["order"] = lineage.get("order", "")
            assembly["family"] = lineage.get("family", "")
            assembly["genus"] = lineage.get("genus", "")

    return assemblies


def write_csv(assemblies: list[dict], output_path: Path):
    """Write assemblies to CSV."""
    if not assemblies:
        print("No assemblies to write!")
        return

    # Sort by taxonomy then N50
    assemblies.sort(key=lambda x: (
        x.get("infraclass", ""),
        x.get("order", ""),
        x.get("family", ""),
        -x.get("scaffold_n50", 0)
    ))

    fieldnames = [
        "accession", "organism_name", "common_name", "species_taxid",
        "infraclass", "order", "family", "genus",
        "assembly_name", "assembly_level", "refseq_category",
        "scaffold_n50", "contig_n50", "total_sequence_length",
        "scaffold_count", "contig_count", "gc_percent", "genome_coverage",
        "gene_count", "protein_coding_count",
        "submission_date", "annotation_name"
    ]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(assemblies)

    print(f"\nWrote {len(assemblies)} assemblies to {output_path}")


def print_summary(assemblies: list[dict]):
    """Print summary statistics."""
    if not assemblies:
        return

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    # Count by major group
    by_infraclass = {}
    by_order = {}

    for a in assemblies:
        ic = a.get("infraclass") or "Unknown"
        order = a.get("order") or "Unknown"
        by_infraclass[ic] = by_infraclass.get(ic, 0) + 1
        by_order[order] = by_order.get(order, 0) + 1

    print(f"\nTotal assemblies: {len(assemblies)}")

    print("\nBy infraclass/class:")
    for ic, count in sorted(by_infraclass.items(), key=lambda x: -x[1])[:10]:
        print(f"  {ic}: {count}")

    print("\nTop orders:")
    for order, count in sorted(by_order.items(), key=lambda x: -x[1])[:15]:
        print(f"  {order}: {count}")

    # N50 stats
    n50s = [a["scaffold_n50"] for a in assemblies]
    print(f"\nScaffold N50 range: {min(n50s)/1e6:.1f} - {max(n50s)/1e6:.1f} Mb")
    print(f"Median N50: {sorted(n50s)[len(n50s)//2]/1e6:.1f} Mb")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch genome assembly quality metrics from NCBI"
    )
    parser.add_argument(
        "--taxon", default="Vertebrata",
        help="Taxonomic group to search (default: Vertebrata)"
    )
    parser.add_argument(
        "--min-n50", type=int, default=50_000_000,
        help="Minimum scaffold N50 in bp (default: 50000000)"
    )
    parser.add_argument(
        "--output", "-o", type=Path,
        default=Path("high_quality_genomes.csv"),
        help="Output CSV path"
    )
    parser.add_argument(
        "--include-genbank", action="store_true",
        help="Include GenBank assemblies (default: RefSeq only)"
    )
    parser.add_argument(
        "--limit", type=int,
        help="Limit number of assemblies (for testing)"
    )
    parser.add_argument(
        "--api-key",
        help="NCBI API key for higher rate limits"
    )
    parser.add_argument(
        "--skip-taxonomy", action="store_true",
        help="Skip taxonomy enrichment (faster)"
    )

    args = parser.parse_args()

    # Fetch assemblies
    assemblies = fetch_assembly_summaries(
        taxon=args.taxon,
        min_n50=args.min_n50,
        refseq_only=not args.include_genbank,
        limit=args.limit,
        api_key=args.api_key
    )

    if not assemblies:
        print("No assemblies found matching criteria!")
        return

    # Enrich with taxonomy
    if not args.skip_taxonomy:
        assemblies = enrich_with_taxonomy(assemblies)

    # Write output
    write_csv(assemblies, args.output)

    # Print summary
    print_summary(assemblies)


if __name__ == "__main__":
    main()
