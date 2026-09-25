#!/usr/bin/env python3
"""
Validate gene copy numbers using miniprot alignment to genome assemblies.
Downloads only the specific genomes needed, runs miniprot, counts loci.

Requirements:
    pip install requests
    brew install miniprot  # or conda install -c bioconda miniprot

Usage:
    python miniprot_validate.py --gene MINDY3 --protein Q9H8M7 \
        --species "Myotis_myotis,Homo_sapiens,Mus_musculus" \
        --output results/miniprot_MINDY3.csv
"""

import argparse
import subprocess
import requests
import gzip
import shutil
import csv
import re
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import time

# NCBI API settings
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
NCBI_DELAY = 0.4


def fetch_protein_sequence(uniprot_id: str) -> Optional[str]:
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        lines = resp.text.strip().split('\n')
        return ''.join(lines[1:])  # Skip header
    except Exception as e:
        print(f"Error fetching protein {uniprot_id}: {e}")
        return None


def get_assembly_accession(species_name: str, session: requests.Session) -> Optional[str]:
    """Find the best genome assembly accession for a species."""
    # Clean species name
    clean_name = species_name.replace("_", " ").split()[0:2]
    search_term = " ".join(clean_name)

    time.sleep(NCBI_DELAY)

    # Search Assembly database
    url = f"{NCBI_BASE}/esearch.fcgi"
    params = {
        "db": "assembly",
        "term": f'"{search_term}"[Organism] AND "latest refseq"[filter]',
        "retmax": 1,
        "retmode": "json"
    }

    try:
        resp = session.get(url, params=params, timeout=30)
        data = resp.json()

        if not data.get("esearchresult", {}).get("idlist"):
            # Try without refseq filter
            params["term"] = f'"{search_term}"[Organism] AND "representative genome"[filter]'
            time.sleep(NCBI_DELAY)
            resp = session.get(url, params=params, timeout=30)
            data = resp.json()

        if not data.get("esearchresult", {}).get("idlist"):
            return None

        assembly_id = data["esearchresult"]["idlist"][0]

        # Get assembly details
        time.sleep(NCBI_DELAY)
        summary_url = f"{NCBI_BASE}/esummary.fcgi"
        summary_params = {
            "db": "assembly",
            "id": assembly_id,
            "retmode": "json"
        }
        resp = session.get(summary_url, params=summary_params, timeout=30)
        summary = resp.json()

        doc = summary.get("result", {}).get(assembly_id, {})
        # Prefer RefSeq, fall back to GenBank
        accession = doc.get("assemblyaccession") or doc.get("gbuid")
        return accession

    except Exception as e:
        print(f"  Error finding assembly for {species_name}: {e}")
        return None


def download_genome(accession: str, output_dir: Path, session: requests.Session) -> Optional[Path]:
    """Download genome assembly using NCBI Datasets."""
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = output_dir / f"{accession}.fna"

    if genome_path.exists():
        print(f"  Using cached genome: {genome_path}")
        return genome_path

    print(f"  Downloading genome {accession}...")

    # Use NCBI Datasets CLI if available
    try:
        # First try datasets CLI
        cmd = [
            "datasets", "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(output_dir / f"{accession}.zip")
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode == 0:
            # Extract the genome
            import zipfile
            zip_path = output_dir / f"{accession}.zip"
            with zipfile.ZipFile(zip_path, 'r') as zf:
                for name in zf.namelist():
                    if name.endswith('.fna') or name.endswith('.fasta'):
                        with zf.open(name) as src, open(genome_path, 'wb') as dst:
                            dst.write(src.read())
                        break
            zip_path.unlink()  # Clean up zip
            return genome_path

    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass

    # Fallback: direct FTP download
    try:
        # Get FTP path from assembly summary
        time.sleep(NCBI_DELAY)
        url = f"{NCBI_BASE}/esearch.fcgi"
        params = {"db": "assembly", "term": accession, "retmode": "json"}
        resp = session.get(url, params=params, timeout=30)
        data = resp.json()

        if data.get("esearchresult", {}).get("idlist"):
            assembly_id = data["esearchresult"]["idlist"][0]

            time.sleep(NCBI_DELAY)
            summary_url = f"{NCBI_BASE}/esummary.fcgi"
            summary_params = {"db": "assembly", "id": assembly_id, "retmode": "json"}
            resp = session.get(summary_url, params=summary_params, timeout=30)
            summary = resp.json()

            doc = summary.get("result", {}).get(assembly_id, {})
            ftp_path = doc.get("ftppath_refseq") or doc.get("ftppath_genbank")

            if ftp_path:
                # Convert FTP to HTTPS
                ftp_path = ftp_path.replace("ftp://", "https://")
                genome_url = f"{ftp_path}/{ftp_path.split('/')[-1]}_genomic.fna.gz"

                print(f"    Downloading from {genome_url}...")
                resp = requests.get(genome_url, stream=True, timeout=600)
                resp.raise_for_status()

                gz_path = output_dir / f"{accession}.fna.gz"
                with open(gz_path, 'wb') as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)

                # Decompress
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(genome_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                gz_path.unlink()

                return genome_path

    except Exception as e:
        print(f"    Download failed: {e}")

    return None


def run_miniprot(protein_fasta: Path, genome_fasta: Path, output_dir: Path) -> Optional[Path]:
    """Run miniprot to align protein to genome."""
    output_gff = output_dir / f"{genome_fasta.stem}_miniprot.gff"

    if output_gff.exists():
        print(f"  Using cached miniprot output: {output_gff}")
        return output_gff

    print(f"  Running miniprot...")

    try:
        cmd = [
            "miniprot",
            "-t", "4",           # threads
            "--gff",             # output GFF
            str(genome_fasta),
            str(protein_fasta)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            print(f"    miniprot error: {result.stderr}")
            return None

        with open(output_gff, 'w') as f:
            f.write(result.stdout)

        return output_gff

    except FileNotFoundError:
        print("    ERROR: miniprot not found. Install with: brew install miniprot")
        return None
    except subprocess.TimeoutExpired:
        print("    ERROR: miniprot timed out")
        return None


def parse_miniprot_gff(gff_path: Path, min_identity: float = 80.0) -> Dict:
    """Parse miniprot GFF output to count gene copies."""
    loci = []

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'mRNA':  # Main alignment record
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                score = float(parts[5]) if parts[5] != '.' else 0
                strand = parts[6]
                attributes = parts[8]

                # Extract identity from attributes
                identity_match = re.search(r'Identity=([0-9.]+)', attributes)
                identity = float(identity_match.group(1)) * 100 if identity_match else 0

                # Extract coverage
                coverage_match = re.search(r'Positive=([0-9.]+)', attributes)
                coverage = float(coverage_match.group(1)) * 100 if coverage_match else 0

                if identity >= min_identity:
                    loci.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'identity': identity,
                        'coverage': coverage,
                        'score': score
                    })

    # Merge overlapping loci (within 50kb on same chromosome/strand)
    merged = []
    loci.sort(key=lambda x: (x['chrom'], x['strand'], x['start']))

    for locus in loci:
        if not merged:
            merged.append(locus)
            continue

        last = merged[-1]
        if (locus['chrom'] == last['chrom'] and
            locus['strand'] == last['strand'] and
            locus['start'] - last['end'] < 50000):
            # Merge: extend end, keep best identity
            last['end'] = max(last['end'], locus['end'])
            last['identity'] = max(last['identity'], locus['identity'])
            last['coverage'] = max(last['coverage'], locus['coverage'])
        else:
            merged.append(locus)

    return {
        'n_loci': len(merged),
        'loci': merged,
        'mean_identity': sum(l['identity'] for l in merged) / len(merged) if merged else 0,
        'best_identity': max((l['identity'] for l in merged), default=0)
    }


def main():
    parser = argparse.ArgumentParser(description="Validate gene copies with miniprot")
    parser.add_argument("--gene", required=True, help="Gene symbol")
    parser.add_argument("--protein", required=True, help="UniProt ID for protein sequence")
    parser.add_argument("--species", required=True, help="Comma-separated species names")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--genome-dir", default="genomes", help="Directory to store genomes")
    parser.add_argument("--min-identity", type=float, default=80.0, help="Minimum identity threshold")
    parser.add_argument("--keep-genomes", action="store_true", help="Keep downloaded genomes")

    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genome_dir = Path(args.genome_dir)

    # Fetch protein sequence
    print(f"\nFetching protein sequence for {args.protein}...")
    protein_seq = fetch_protein_sequence(args.protein)
    if not protein_seq:
        print("ERROR: Could not fetch protein sequence")
        return 1

    protein_fasta = genome_dir / f"{args.protein}.fasta"
    genome_dir.mkdir(parents=True, exist_ok=True)
    with open(protein_fasta, 'w') as f:
        f.write(f">{args.protein}\n{protein_seq}\n")
    print(f"  Protein length: {len(protein_seq)} aa")

    # Process each species
    species_list = [s.strip() for s in args.species.split(',')]
    results = []

    session = requests.Session()

    for species in species_list:
        print(f"\n{'='*60}")
        print(f"Processing: {species}")
        print('='*60)

        result = {
            'gene': args.gene,
            'species': species,
            'status': 'failed',
            'assembly': '',
            'n_loci': 0,
            'mean_identity': 0,
            'best_identity': 0,
            'loci_details': ''
        }

        # Find assembly
        accession = get_assembly_accession(species, session)
        if not accession:
            print(f"  No assembly found for {species}")
            results.append(result)
            continue

        result['assembly'] = accession
        print(f"  Assembly: {accession}")

        # Download genome
        genome_path = download_genome(accession, genome_dir, session)
        if not genome_path:
            print(f"  Failed to download genome")
            results.append(result)
            continue

        # Run miniprot
        gff_path = run_miniprot(protein_fasta, genome_path, genome_dir)
        if not gff_path:
            results.append(result)
            continue

        # Parse results
        parse_result = parse_miniprot_gff(gff_path, args.min_identity)

        result['status'] = 'success'
        result['n_loci'] = parse_result['n_loci']
        result['mean_identity'] = round(parse_result['mean_identity'], 2)
        result['best_identity'] = round(parse_result['best_identity'], 2)
        result['loci_details'] = ';'.join(
            f"{l['chrom']}:{l['start']}-{l['end']}({l['identity']:.1f}%)"
            for l in parse_result['loci'][:5]  # Top 5 loci
        )

        print(f"  Found {result['n_loci']} loci (identity >= {args.min_identity}%)")
        for locus in parse_result['loci'][:3]:
            print(f"    - {locus['chrom']}:{locus['start']}-{locus['end']} "
                  f"({locus['strand']}, {locus['identity']:.1f}% identity)")

        results.append(result)

        # Clean up genome if not keeping
        if not args.keep_genomes and genome_path.exists():
            genome_path.unlink()
            if gff_path.exists():
                gff_path.unlink()

    # Write results
    print(f"\n{'='*60}")
    print(f"Writing results to {output_path}")
    print('='*60)

    fieldnames = ['gene', 'species', 'status', 'assembly', 'n_loci',
                  'mean_identity', 'best_identity', 'loci_details']

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # Summary
    print("\nSummary:")
    print("-" * 40)
    for r in results:
        status = "✓" if r['status'] == 'success' else "✗"
        print(f"  {status} {r['species']}: {r['n_loci']} loci")

    return 0


if __name__ == "__main__":
    exit(main())
