#!/usr/bin/env python3
"""
Per-gene miniprot validation for TOGA loss/duplication candidates.

Creates individual output files for each gene: results/{gene}_validation.csv

Key features:
1. Skips genes where no affected species has an available genome
2. Includes affected species + outgroup controls for comparison
3. Outputs per-gene files for easy searching

Usage:
    python validate_per_gene.py \
        --candidates results/TOP_CANDIDATES/refined_candidates_for_miniprot.csv \
        --toga-loss-matrix data/loss_analysis/All_Species_Gene_PresenceAbsence.tsv \
        --toga-dup-matrix data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv \
        --output-dir results/gene_validations \
        --max-genes 10
"""

import argparse
import subprocess
import requests
import gzip
import shutil
import csv
import sys
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Set
import time

# NCBI API settings
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_DELAY = 0.4

# Reference species (always include for baseline comparison)
REFERENCE_SPECIES = ["Homo_sapiens", "Mus_musculus"]


def fetch_protein_sequence(gene_symbol: str, session: requests.Session) -> Tuple[Optional[str], Optional[str]]:
    """Fetch canonical protein sequence from UniProt by gene symbol."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    # Prefer reviewed (Swiss-Prot) entries for canonical sequences
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "fasta",
        "size": 1
    }
    try:
        resp = session.get(url, params=params, timeout=30)
        resp.raise_for_status()
        if resp.text.strip():
            lines = resp.text.strip().split('\n')
            header = lines[0]
            uniprot_id = header.split('|')[1] if '|' in header else None
            sequence = ''.join(lines[1:])
            return uniprot_id, sequence
        # Fall back to unreviewed if no reviewed entry
        params["query"] = f"gene:{gene_symbol} AND organism_id:9606"
        resp = session.get(url, params=params, timeout=30)
        resp.raise_for_status()
        if resp.text.strip():
            lines = resp.text.strip().split('\n')
            header = lines[0]
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
        accession = doc.get("assemblyaccession", "")
        ftp_path = doc.get("ftppath_refseq") or doc.get("ftppath_genbank", "")

        try:
            n50 = int(doc.get("scaffoldn50", 0) or 0)
        except:
            n50 = 0

        return {
            "accession": accession,
            "scaffold_n50": n50,
            "ftp_path": ftp_path
        }
    except Exception as e:
        print(f"    Error getting assembly for {species_name}: {e}")
        return None


def download_genome(accession: str, ftp_path: str, genome_dir: Path) -> Optional[Path]:
    """Download genome FASTA from NCBI."""
    local_path = genome_dir / f"{accession}.fna"
    if local_path.exists():
        return local_path

    if not ftp_path:
        return None

    filename = ftp_path.split("/")[-1]
    # Convert FTP to HTTPS URL (requests can't handle FTP)
    https_path = ftp_path.replace("ftp://", "https://")
    fasta_url = f"{https_path}/{filename}_genomic.fna.gz"
    gz_path = genome_dir / f"{accession}.fna.gz"

    try:
        print(f"    Downloading {accession}...", flush=True)
        resp = requests.get(fasta_url, stream=True, timeout=300)
        resp.raise_for_status()
        with open(gz_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)

        with gzip.open(gz_path, 'rb') as f_in:
            with open(local_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        gz_path.unlink()
        return local_path
    except Exception as e:
        print(f"    Error downloading {accession}: {e}")
        if gz_path.exists():
            gz_path.unlink()
        return None


def run_miniprot(protein_fasta: Path, genome_fasta: Path, output_gff: Path) -> bool:
    """Run miniprot to align protein to genome."""
    try:
        cmd = ["miniprot", "-t4", "--gff", str(genome_fasta), str(protein_fasta)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=2700)  # 45 min for large genomes
        with open(output_gff, 'w') as f:
            f.write(result.stdout)
        return True
    except subprocess.TimeoutExpired:
        print(f"    miniprot timeout (45 min) on {genome_fasta.name}", flush=True)
        return False
    except Exception as e:
        print(f"    miniprot error: {e}")
        return False


def parse_miniprot_gff(gff_path: Path, min_identity: float = 50.0) -> Dict:
    """Parse miniprot GFF output to count distinct loci."""
    loci = []

    if not gff_path.exists():
        return {'n_loci': 0, 'n_raw_loci': 0, 'best_identity': 0, 'loci': []}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue

            attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
            identity = float(attrs.get('Identity', 0))

            if identity >= min_identity / 100:  # Identity is 0-1 in GFF
                loci.append({
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'identity': identity * 100  # Convert to percentage
                })

    n_raw_loci = len(loci)

    # Merge overlapping loci on same chromosome/strand
    loci.sort(key=lambda x: (x['chrom'], x['strand'], x['start']))
    merged = []
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
        'n_raw_loci': n_raw_loci,
        'best_identity': max((l['identity'] for l in merged), default=0),
        'loci': merged
    }


def load_toga_matrix(matrix_path: Path) -> Dict[str, Dict[str, int]]:
    """Load TOGA binary matrix (presence/absence or duplication)."""
    toga_data = {}
    with open(matrix_path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        species_cols = header[2:]  # Skip gene_id, gene_symbol

        for row in reader:
            if len(row) < 3:
                continue
            gene_symbol = row[1]
            toga_data[gene_symbol] = {}
            for i, species in enumerate(species_cols):
                species_clean = species.replace(".", "_")
                try:
                    toga_data[gene_symbol][species_clean] = int(row[i + 2])
                except (ValueError, IndexError):
                    pass  # Leave as missing

    return toga_data


def load_refined_candidates(candidates_path: Path) -> List[Dict]:
    """Load refined candidate genes."""
    candidates = []
    with open(candidates_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append({
                'gene': row.get('gene', ''),
                'gene_id': row.get('gene_id', ''),
                'event_type': row.get('event_type', 'unknown'),
                'source': row.get('source', 'unknown'),
                'priority': int(row.get('priority', 99)),
                'priority_reason': row.get('priority_reason', ''),
                'n_events': int(row.get('n_events', 0)) if row.get('n_events') else 0,
                'n_clades': int(row.get('n_clades', 0)) if row.get('n_clades') else 0,
                'is_longevity_gene': row.get('is_longevity_gene', '').lower() == 'true'
            })
    # Sort by priority
    candidates.sort(key=lambda x: x['priority'])
    return candidates


def get_affected_and_control_species(gene: str, event_type: str,
                                      loss_matrix: Dict, dup_matrix: Dict) -> Tuple[List[str], List[str]]:
    """Get species with and without the event from TOGA matrix."""
    if event_type == 'loss':
        gene_data = loss_matrix.get(gene, {})
        affected = [sp for sp, val in gene_data.items() if val == 0]
        controls = [sp for sp, val in gene_data.items() if val == 1]
    else:  # duplication
        gene_data = dup_matrix.get(gene, {})
        affected = [sp for sp, val in gene_data.items() if val == 1]
        controls = [sp for sp, val in gene_data.items() if val == 0]

    return affected, controls


def main():
    parser = argparse.ArgumentParser(description="Per-gene miniprot validation")
    parser.add_argument("--candidates", required=True,
                        help="CSV with candidate genes (refined_candidates_for_miniprot.csv)")
    parser.add_argument("--toga-loss-matrix", required=True,
                        help="TOGA presence/absence matrix (0=lost, 1=present)")
    parser.add_argument("--toga-dup-matrix", required=True,
                        help="TOGA duplication binary matrix (1=duplicated)")
    parser.add_argument("--output-dir", default="results/gene_validations",
                        help="Directory for per-gene output files")
    parser.add_argument("--genome-dir", default="genomes", help="Directory for genome files")
    parser.add_argument("--min-n50", type=int, default=1000000,
                        help="Minimum scaffold N50 for quality filter")
    parser.add_argument("--min-identity", type=float, default=50.0,
                        help="Minimum protein identity")
    parser.add_argument("--max-genes", type=int, default=20,
                        help="Maximum genes to validate")
    parser.add_argument("--max-species-per-gene", type=int, default=8,
                        help="Maximum species per gene")
    parser.add_argument("--keep-genomes", action="store_true",
                        help="Keep downloaded genomes")
    parser.add_argument("--skip-download", action="store_true",
                        help="Only use already-cached genomes")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    genome_dir = Path(args.genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    # Load candidates
    print("Loading candidates...", flush=True)
    candidates = load_refined_candidates(Path(args.candidates))
    print(f"  Loaded {len(candidates)} candidates (sorted by priority)", flush=True)

    # Load TOGA matrices
    print("\nLoading TOGA matrices...", flush=True)
    loss_matrix = load_toga_matrix(Path(args.toga_loss_matrix))
    print(f"  Loss matrix: {len(loss_matrix)} genes", flush=True)

    dup_matrix = load_toga_matrix(Path(args.toga_dup_matrix))
    print(f"  Duplication matrix: {len(dup_matrix)} genes", flush=True)

    # Pre-check genome availability for ALL candidate species (not just first N)
    print("\nChecking genome assemblies for all species...", flush=True)
    all_species = set()
    for cand in candidates:  # Check ALL candidates
        affected, controls = get_affected_and_control_species(
            cand['gene'], cand['event_type'], loss_matrix, dup_matrix
        )
        all_species.update(affected[:4])  # Up to 4 affected
        all_species.update(REFERENCE_SPECIES)
        all_species.update(controls[:4])  # Up to 4 controls

    species_info = {}
    genome_paths = {}

    for species in sorted(all_species):
        info = get_assembly_info(species, session)
        if info and info.get('scaffold_n50', 0) >= args.min_n50:
            species_info[species] = info
            # Check if cached
            cached_path = genome_dir / f"{info['accession']}.fna"
            if cached_path.exists():
                genome_paths[species] = cached_path
                print(f"  ✓ {species}: N50={info['scaffold_n50']:,} [cached]", flush=True)
            else:
                print(f"  ✓ {species}: N50={info['scaffold_n50']:,}", flush=True)
        elif info:
            print(f"  ✗ {species}: N50={info.get('scaffold_n50', 0):,} (below threshold)", flush=True)
        else:
            print(f"  ✗ {species}: No assembly found", flush=True)

    # Filter ALL candidates to those with at least one affected species genome available
    # Then limit to max_genes
    print("\nFiltering candidates with available affected species genomes...", flush=True)
    valid_candidates = []
    skipped = 0
    for cand in candidates:  # Check ALL candidates
        affected, _ = get_affected_and_control_species(
            cand['gene'], cand['event_type'], loss_matrix, dup_matrix
        )
        affected_with_genome = [sp for sp in affected if sp in species_info]
        if affected_with_genome:
            cand['affected_with_genome'] = affected_with_genome
            valid_candidates.append(cand)
            if len(valid_candidates) <= args.max_genes:
                print(f"  ✓ {cand['gene']}: {len(affected_with_genome)} affected species with genomes", flush=True)
        else:
            skipped += 1

    print(f"  ... {skipped} candidates skipped (no affected species genomes)", flush=True)

    # Limit to max_genes
    valid_candidates = valid_candidates[:args.max_genes]

    print(f"\n{len(valid_candidates)} genes with validatable affected species", flush=True)

    if not valid_candidates:
        print("ERROR: No genes can be validated (no affected species have available genomes)")
        return 1

    # Download needed genomes
    if not args.skip_download:
        print("\nDownloading genomes...", flush=True)
        for species, info in species_info.items():
            if species not in genome_paths and info['ftp_path']:
                path = download_genome(info['accession'], info['ftp_path'], genome_dir)
                if path:
                    genome_paths[species] = path

    # Validate each gene
    for cand in valid_candidates:
        gene = cand['gene']
        event_type = cand['event_type']

        print(f"\n{'='*70}", flush=True)
        print(f"Gene: {gene} ({event_type.upper()})", flush=True)
        print(f"Priority: {cand['priority']} - {cand['priority_reason']}", flush=True)
        print(f"Events: {cand['n_events']} in {cand['n_clades']} clades", flush=True)
        print('='*70, flush=True)

        # Get protein sequence
        uniprot_id, protein_seq = fetch_protein_sequence(gene, session)
        if not protein_seq:
            print(f"  Could not fetch protein sequence - skipping", flush=True)
            continue

        protein_fasta = genome_dir / f"{gene}.fasta"
        with open(protein_fasta, 'w') as f:
            f.write(f">{gene}\n{protein_seq}\n")
        print(f"  Protein: {uniprot_id}, {len(protein_seq)} aa", flush=True)

        # Get affected and control species
        affected, controls = get_affected_and_control_species(gene, event_type, loss_matrix, dup_matrix)

        # Select species to validate
        species_to_check = []

        # Add affected species with genomes (up to 4)
        for sp in affected:
            if sp in genome_paths and sp not in species_to_check:
                species_to_check.append(sp)
                if len([s for s in species_to_check if s in affected]) >= 4:
                    break

        # Add reference species as controls
        for ref in REFERENCE_SPECIES:
            if ref in genome_paths and ref not in species_to_check:
                species_to_check.append(ref)

        # Add more controls if needed
        for sp in controls:
            if sp in genome_paths and sp not in species_to_check:
                species_to_check.append(sp)
                if len(species_to_check) >= args.max_species_per_gene:
                    break

        print(f"  Validating {len([s for s in species_to_check if s in affected])} affected + "
              f"{len([s for s in species_to_check if s in controls or s in REFERENCE_SPECIES])} controls", flush=True)

        # Get the appropriate TOGA matrix
        if event_type == 'loss':
            toga_data = loss_matrix
            toga_field = 'toga_presence'
        else:
            toga_data = dup_matrix
            toga_field = 'toga_duplication'

        # Per-gene output file
        gene_output = output_dir / f"{gene}_validation.csv"
        fieldnames = ['gene', 'event_type', 'priority', 'priority_reason', 'n_events', 'n_clades',
                      'species', 'is_affected_species', 'assembly', 'scaffold_n50',
                      'miniprot_loci', 'miniprot_raw_loci', 'miniprot_best_identity',
                      toga_field, 'concordant', 'validation_result']

        results = []

        with open(gene_output, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()
            outfile.flush()

            for species in species_to_check:
                if species not in genome_paths:
                    print(f"    ✗ {species}: genome not available", flush=True)
                    continue

                genome_path = genome_paths[species]
                gff_path = genome_dir / f"{species}_{gene}.gff"

                # Run miniprot
                print(f"    Running miniprot on {species}...", flush=True)
                if not run_miniprot(protein_fasta, genome_path, gff_path):
                    continue

                # Parse results
                mp_result = parse_miniprot_gff(gff_path, args.min_identity)

                # Get TOGA prediction
                toga_value = toga_data.get(gene, {}).get(species, None)
                if toga_value is None:
                    # Try partial match
                    gene_data = toga_data.get(gene, {})
                    for toga_species, value in gene_data.items():
                        if toga_species.startswith(species) or species.startswith(toga_species.rsplit('_', 1)[0]):
                            toga_value = value
                            break

                # Determine concordance and validation result
                is_affected = species in affected

                if event_type == 'loss':
                    # Loss validation
                    if is_affected:
                        # Affected species: TOGA says lost (0), expect 0 loci
                        concordant = 'yes' if mp_result['n_loci'] == 0 else 'no'
                        validation = 'TRUE_LOSS' if mp_result['n_loci'] == 0 else 'FALSE_LOSS'
                    else:
                        # Control species: TOGA says present (1), expect >= 1 loci
                        concordant = 'yes' if mp_result['n_loci'] > 0 else 'no'
                        validation = 'PRESENT' if mp_result['n_loci'] > 0 else 'MISSING_IN_CONTROL'
                else:
                    # Duplication validation
                    if is_affected:
                        # Affected species: TOGA says duplicated (1), expect > 1 loci
                        concordant = 'yes' if mp_result['n_loci'] > 1 else 'no'
                        validation = 'TRUE_DUP' if mp_result['n_loci'] > 1 else 'FALSE_DUP'
                    else:
                        # Control species: TOGA says single copy (0), expect 1 locus
                        concordant = 'yes' if mp_result['n_loci'] == 1 else 'no'
                        validation = 'SINGLE_COPY' if mp_result['n_loci'] == 1 else 'UNEXPECTED'

                result = {
                    'gene': gene,
                    'event_type': event_type,
                    'priority': cand['priority'],
                    'priority_reason': cand['priority_reason'],
                    'n_events': cand['n_events'],
                    'n_clades': cand['n_clades'],
                    'species': species,
                    'is_affected_species': 'yes' if is_affected else 'no',
                    'assembly': species_info.get(species, {}).get('accession', ''),
                    'scaffold_n50': species_info.get(species, {}).get('scaffold_n50', 0),
                    'miniprot_loci': mp_result['n_loci'],
                    'miniprot_raw_loci': mp_result['n_raw_loci'],
                    'miniprot_best_identity': round(mp_result['best_identity'], 1),
                    toga_field: toga_value,
                    'concordant': concordant,
                    'validation_result': validation
                }
                results.append(result)

                # Write result immediately
                writer.writerow(result)
                outfile.flush()

                status = "✓" if concordant == 'yes' else "✗"
                affected_marker = "[AFFECTED]" if is_affected else "[CONTROL]"
                print(f"    {status} {species} {affected_marker}: {mp_result['n_loci']} loci, "
                      f"TOGA={toga_value} -> {validation}", flush=True)

                # Clean up GFF
                if not args.keep_genomes and gff_path.exists():
                    gff_path.unlink()

        # Gene summary
        if results:
            affected_results = [r for r in results if r['is_affected_species'] == 'yes']
            control_results = [r for r in results if r['is_affected_species'] == 'no']

            affected_concordant = sum(1 for r in affected_results if r['concordant'] == 'yes')
            control_concordant = sum(1 for r in control_results if r['concordant'] == 'yes')

            print(f"\n  Summary for {gene}:", flush=True)
            print(f"    Affected species: {affected_concordant}/{len(affected_results)} concordant", flush=True)
            print(f"    Control species: {control_concordant}/{len(control_results)} concordant", flush=True)
            print(f"    Output: {gene_output}", flush=True)

    print(f"\n{'='*70}", flush=True)
    print(f"Validation complete! Per-gene files in: {output_dir}", flush=True)
    print('='*70, flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
