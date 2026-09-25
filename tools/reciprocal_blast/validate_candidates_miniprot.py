#!/usr/bin/env python3
"""
Systematic miniprot validation of TOGA duplication/loss candidates.

Replicates the CGAS/armadillo validation approach:
1. Load candidate genes from analysis results
2. Filter species by genome quality (N50 threshold)
3. Run miniprot to count true genomic loci
4. Compare to TOGA predictions

Usage:
    python validate_candidates_miniprot.py \
        --candidates results/TOP_CANDIDATES/candidates_for_blat.csv \
        --toga-matrix data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv \
        --species-traits data/species_trait_data.csv \
        --output results/miniprot_validation_summary.csv \
        --min-n50 1000000 \
        --max-genes 20
"""

import argparse
import subprocess
import requests
import gzip
import shutil
import csv
import re
import json
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import time

# NCBI API settings
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_DELAY = 0.4


def fetch_protein_sequence(gene_symbol: str, session: requests.Session) -> Tuple[Optional[str], Optional[str]]:
    """Fetch protein sequence from UniProt by gene symbol."""
    # Search UniProt for human protein
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606",
        "format": "fasta",
        "size": 1
    }
    try:
        resp = session.get(url, params=params, timeout=30)
        resp.raise_for_status()
        if resp.text.strip():
            lines = resp.text.strip().split('\n')
            header = lines[0]
            # Extract UniProt ID from header
            uniprot_id = header.split('|')[1] if '|' in header else None
            sequence = ''.join(lines[1:])
            return uniprot_id, sequence
    except Exception as e:
        print(f"    Error fetching protein for {gene_symbol}: {e}")
    return None, None


def get_assembly_info(species_name: str, session: requests.Session) -> Optional[Dict]:
    """Get assembly accession and quality metrics for a species."""
    clean_name = species_name.replace("_", " ").split()[0:2]
    search_term = " ".join(clean_name)

    time.sleep(NCBI_DELAY)

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

        time.sleep(NCBI_DELAY)
        summary_url = f"{NCBI_BASE}/esummary.fcgi"
        summary_params = {"db": "assembly", "id": assembly_id, "retmode": "json"}
        resp = session.get(summary_url, params=summary_params, timeout=30)
        summary = resp.json()

        doc = summary.get("result", {}).get(assembly_id, {})

        return {
            "accession": doc.get("assemblyaccession") or doc.get("gbuid"),
            "assembly_level": doc.get("assemblylevel", ""),
            "contig_n50": int(doc.get("contign50", 0)),
            "scaffold_n50": int(doc.get("scaffoldn50", 0)),
            "ftp_path": doc.get("ftppath_refseq") or doc.get("ftppath_genbank")
        }

    except Exception as e:
        print(f"    Error finding assembly for {species_name}: {e}")
        return None


def download_genome(accession: str, ftp_path: str, output_dir: Path) -> Optional[Path]:
    """Download genome assembly."""
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_path = output_dir / f"{accession}.fna"

    if genome_path.exists():
        return genome_path

    print(f"    Downloading {accession}...")

    try:
        # Convert FTP to HTTPS
        ftp_path = ftp_path.replace("ftp://", "https://")
        genome_url = f"{ftp_path}/{ftp_path.split('/')[-1]}_genomic.fna.gz"

        resp = requests.get(genome_url, stream=True, timeout=600)
        resp.raise_for_status()

        gz_path = output_dir / f"{accession}.fna.gz"
        with open(gz_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)

        with gzip.open(gz_path, 'rb') as f_in:
            with open(genome_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        gz_path.unlink()

        return genome_path

    except Exception as e:
        print(f"    Download failed: {e}")
        return None


def run_miniprot(protein_fasta: Path, genome_fasta: Path, output_gff: Path) -> bool:
    """Run miniprot alignment."""
    if output_gff.exists():
        return True

    try:
        cmd = ["miniprot", "-t", "4", "--gff", str(genome_fasta), str(protein_fasta)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            return False

        with open(output_gff, 'w') as f:
            f.write(result.stdout)
        return True

    except Exception as e:
        print(f"    miniprot error: {e}")
        return False


def parse_miniprot_gff(gff_path: Path, min_identity: float = 50.0) -> Dict:
    """Parse miniprot GFF to count distinct genomic loci."""
    loci = []

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            identity_match = re.search(r'Identity=([0-9.]+)', attributes)
            identity = float(identity_match.group(1)) * 100 if identity_match else 0

            if identity >= min_identity:
                loci.append({
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'identity': identity
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
            last['end'] = max(last['end'], locus['end'])
            last['identity'] = max(last['identity'], locus['identity'])
        else:
            merged.append(locus)

    return {
        'n_loci': len(merged),
        'best_identity': max((l['identity'] for l in merged), default=0),
        'loci': merged
    }


def load_toga_matrix(matrix_path: Path) -> Dict[str, Dict[str, int]]:
    """Load TOGA duplication binary matrix."""
    toga_data = {}
    with open(matrix_path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        species_cols = header[2:]  # Skip gene_id, gene_symbol

        for row in reader:
            gene_id = row[0]
            gene_symbol = row[1]
            toga_data[gene_symbol] = {}
            for i, species in enumerate(species_cols):
                species_clean = species.replace(".", "_")
                try:
                    toga_data[gene_symbol][species_clean] = int(row[i + 2])
                except (ValueError, IndexError):
                    toga_data[gene_symbol][species_clean] = 0

    return toga_data


def main():
    parser = argparse.ArgumentParser(description="Validate TOGA candidates with miniprot")
    parser.add_argument("--candidates", help="CSV with candidate genes (gene_symbol column required)")
    parser.add_argument("--genes", help="Comma-separated gene symbols (alternative to --candidates)")
    parser.add_argument("--toga-matrix", required=True, help="TOGA presence/absence matrix for loss validation")
    parser.add_argument("--species", help="Comma-separated species to validate")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--genome-dir", default="genomes", help="Directory for genome files")
    parser.add_argument("--min-n50", type=int, default=1000000, help="Minimum scaffold N50 for quality filter")
    parser.add_argument("--min-identity", type=float, default=50.0, help="Minimum protein identity")
    parser.add_argument("--max-genes", type=int, default=20, help="Maximum genes to validate")
    parser.add_argument("--keep-genomes", action="store_true", help="Keep downloaded genomes")
    parser.add_argument("--skip-download", action="store_true", help="Only use already-cached genomes")

    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genome_dir = Path(args.genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    # Load candidate genes
    if args.genes:
        genes = [g.strip() for g in args.genes.split(',')]
    elif args.candidates:
        genes = []
        with open(args.candidates) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'gene_symbol' in row:
                    genes.append(row['gene_symbol'])
        genes = genes[:args.max_genes]
    else:
        print("ERROR: Must provide --candidates or --genes")
        return 1

    print(f"Validating {len(genes)} genes: {', '.join(genes[:5])}...")

    # Load TOGA matrix
    print(f"\nLoading TOGA matrix from {args.toga_matrix}...")
    toga_data = load_toga_matrix(Path(args.toga_matrix))
    print(f"  Loaded {len(toga_data)} genes")

    # Parse species list
    if args.species:
        species_list = [s.strip() for s in args.species.split(',')]
    else:
        # Default: use species from TOGA matrix for first gene
        first_gene = genes[0] if genes else None
        if first_gene and first_gene in toga_data:
            species_list = list(toga_data[first_gene].keys())[:10]
        else:
            species_list = ["Homo_sapiens", "Mus_musculus", "Myotis_myotis"]

    print(f"Species to validate: {len(species_list)}")

    # Get assembly info for all species (with quality filtering)
    print("\nChecking genome assemblies...")
    species_info = {}
    for species in species_list:
        info = get_assembly_info(species, session)
        if info and info.get('scaffold_n50', 0) >= args.min_n50:
            species_info[species] = info
            print(f"  ✓ {species}: N50={info['scaffold_n50']:,}")
        elif info:
            print(f"  ✗ {species}: N50={info.get('scaffold_n50', 0):,} (below threshold)")
        else:
            print(f"  ✗ {species}: No assembly found")

    if not species_info:
        print("ERROR: No species passed quality filter")
        return 1

    # Download genomes (or use cached only)
    if args.skip_download:
        print("\nUsing cached genomes only (--skip-download)...")
        genome_paths = {}
        # Find all cached .fna files
        cached_files = list(genome_dir.glob("*.fna"))
        cached_accessions = {f.stem: f for f in cached_files}
        for species, info in species_info.items():
            acc = info['accession']
            if acc in cached_accessions:
                genome_paths[species] = cached_accessions[acc]
                print(f"  ✓ {species}: using cached {acc}")
            else:
                print(f"  ✗ {species}: not cached (skipping)")
    else:
        print("\nDownloading genomes...")
        genome_paths = {}
        for species, info in species_info.items():
            if info['ftp_path']:
                path = download_genome(info['accession'], info['ftp_path'], genome_dir)
                if path:
                    genome_paths[species] = path

    # Validate each gene
    results = []

    for gene in genes:
        print(f"\n{'='*60}")
        print(f"Gene: {gene}")
        print('='*60)

        # Get protein sequence
        uniprot_id, protein_seq = fetch_protein_sequence(gene, session)
        if not protein_seq:
            print(f"  Could not fetch protein sequence")
            continue

        protein_fasta = genome_dir / f"{gene}.fasta"
        with open(protein_fasta, 'w') as f:
            f.write(f">{gene}\n{protein_seq}\n")
        print(f"  Protein: {uniprot_id}, {len(protein_seq)} aa")

        # Validate in each species
        for species, genome_path in genome_paths.items():
            gff_path = genome_dir / f"{species}_{gene}.gff"

            # Run miniprot
            if not run_miniprot(protein_fasta, genome_path, gff_path):
                continue

            # Parse results
            mp_result = parse_miniprot_gff(gff_path, args.min_identity)

            # Get TOGA prediction (presence/absence: 1=present, 0=lost)
            # Handle species name suffixes (e.g., Bradypus_torquatus vs Bradypus_torquatus_1)
            toga_presence = toga_data.get(gene, {}).get(species, None)
            if toga_presence is None:
                # Try partial match for species with suffixes
                gene_data = toga_data.get(gene, {})
                for toga_species, value in gene_data.items():
                    if toga_species.startswith(species) or species.startswith(toga_species.rsplit('_', 1)[0]):
                        toga_presence = value
                        break

            # Record result - concordance for LOSS validation:
            # miniprot loci > 0 (gene found) should match TOGA presence == 1
            result = {
                'gene': gene,
                'species': species,
                'assembly': species_info[species]['accession'],
                'scaffold_n50': species_info[species]['scaffold_n50'],
                'miniprot_loci': mp_result['n_loci'],
                'miniprot_best_identity': round(mp_result['best_identity'], 1),
                'toga_presence': toga_presence,
                'concordant': 'yes' if (mp_result['n_loci'] > 0) == (toga_presence == 1) else 'no'
            }
            results.append(result)

            status = "✓" if result['concordant'] == 'yes' else "✗"
            print(f"  {status} {species}: miniprot={mp_result['n_loci']} loci, TOGA_presence={toga_presence}")

            # Clean up GFF if not keeping
            if not args.keep_genomes and gff_path.exists():
                gff_path.unlink()

    # Write results
    print(f"\n{'='*60}")
    print(f"Writing results to {output_path}")
    print('='*60)

    fieldnames = ['gene', 'species', 'assembly', 'scaffold_n50',
                  'miniprot_loci', 'miniprot_best_identity',
                  'toga_presence', 'concordant']

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # Summary
    concordant = sum(1 for r in results if r['concordant'] == 'yes')
    total = len(results)
    print(f"\nConcordance: {concordant}/{total} ({100*concordant/total:.1f}%)")

    return 0


if __name__ == "__main__":
    exit(main())
