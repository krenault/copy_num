#!/usr/bin/env python3
"""
Unified miniprot validation for TOGA loss/duplication candidates.

Key improvements over previous scripts:
1. Includes WHY each gene was flagged (pattern, affected species) in output
2. Auto-selects species with variation + reference controls
3. Ensures meaningful validation by requiring species showing the event

Usage:
    python validate_with_context.py \
        --candidates results/TOP_CANDIDATES/interesting_longevity_candidates.csv \
        --toga-loss-matrix data/loss_analysis/All_Species_Gene_PresenceAbsence.tsv \
        --toga-dup-matrix data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv \
        --output results/validation_with_context.csv \
        --max-genes 10
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
        print(f"    Downloading {accession}...")
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

            if identity >= min_identity:
                loci.append({
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'strand': parts[6],
                    'identity': identity
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


def load_candidates(candidates_path: Path) -> List[Dict]:
    """Load candidate genes with their context (why flagged)."""
    candidates = []
    with open(candidates_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append({
                'gene': row.get('gene', row.get('gene_symbol', '')),
                'gene_id': row.get('gene_id', ''),
                'event_type': row.get('event_type', 'unknown'),
                'pattern': row.get('pattern', 'unknown'),
                'n_events': int(row.get('n_events', 0)),
                'primary_clade': row.get('primary_clade', ''),
                'clade_concentration': float(row.get('clade_concentration', 0)),
                'species_list': row.get('species_list', '').split(';') if row.get('species_list') else []
            })
    return candidates


def get_species_for_gene(gene: str, candidate_info: Dict,
                         loss_matrix: Dict, dup_matrix: Dict,
                         max_species: int = 8) -> Tuple[List[str], str]:
    """
    Select species for validation ensuring variation is included.

    Returns:
        - List of species to validate
        - Reason string explaining selection
    """
    event_type = candidate_info.get('event_type', 'unknown')
    pattern = candidate_info.get('pattern', 'unknown')
    affected_species = candidate_info.get('species_list', [])

    # Get species showing the event from the TOGA matrix
    if event_type == 'loss':
        gene_data = loss_matrix.get(gene, {})
        # For loss: presence/absence matrix where 0 = lost
        species_with_event = [sp for sp, val in gene_data.items() if val == 0]
        species_without_event = [sp for sp, val in gene_data.items() if val == 1]
    else:  # duplication
        gene_data = dup_matrix.get(gene, {})
        # For duplication: 1 = duplicated
        species_with_event = [sp for sp, val in gene_data.items() if val == 1]
        species_without_event = [sp for sp, val in gene_data.items() if val == 0]

    # Use affected_species from candidate info if available
    if affected_species:
        species_with_event = list(set(species_with_event) | set(affected_species))

    # Build selection
    selected = []

    # 1. Always include reference species that don't have the event (for baseline)
    for ref in REFERENCE_SPECIES:
        if ref in species_without_event or ref not in species_with_event:
            if ref not in selected:
                selected.append(ref)

    # 2. Include species WITH the event (this is critical!)
    n_event_species = min(4, len(species_with_event))  # Up to 4 affected species
    for sp in species_with_event[:n_event_species]:
        if sp not in selected:
            selected.append(sp)

    # 3. Add a few more species without the event for comparison
    remaining_slots = max_species - len(selected)
    for sp in species_without_event:
        if sp not in selected and remaining_slots > 0:
            selected.append(sp)
            remaining_slots -= 1

    # Build reason string
    n_affected = len([s for s in selected if s in species_with_event])
    n_controls = len([s for s in selected if s in species_without_event])

    reason = f"{pattern}: {candidate_info.get('n_events', 0)} {event_type}(s)"
    if candidate_info.get('primary_clade'):
        reason += f" in {candidate_info['primary_clade']}"
    reason += f" | validating {n_affected} affected + {n_controls} controls"

    return selected, reason


def main():
    parser = argparse.ArgumentParser(description="Validate TOGA candidates with context")
    parser.add_argument("--candidates", required=True,
                        help="CSV with candidate genes (from find_interesting_candidates.py)")
    parser.add_argument("--toga-loss-matrix", required=True,
                        help="TOGA presence/absence matrix (0=lost, 1=present)")
    parser.add_argument("--toga-dup-matrix", required=True,
                        help="TOGA duplication binary matrix (1=duplicated)")
    parser.add_argument("--output", required=True, help="Output CSV path")
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

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genome_dir = Path(args.genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    # Load candidates with context
    print("Loading candidates...")
    candidates = load_candidates(Path(args.candidates))
    print(f"  Loaded {len(candidates)} candidates")

    # Limit to max_genes
    candidates = candidates[:args.max_genes]
    print(f"  Processing {len(candidates)} genes")

    # Load TOGA matrices
    print("\nLoading TOGA matrices...")
    loss_matrix = load_toga_matrix(Path(args.toga_loss_matrix))
    print(f"  Loss matrix: {len(loss_matrix)} genes")

    dup_matrix = load_toga_matrix(Path(args.toga_dup_matrix))
    print(f"  Duplication matrix: {len(dup_matrix)} genes")

    # Collect all species we need
    all_needed_species: Set[str] = set()
    gene_species_map = {}
    gene_reasons = {}

    print("\nSelecting species for each gene...")
    for cand in candidates:
        gene = cand['gene']
        species_list, reason = get_species_for_gene(
            gene, cand, loss_matrix, dup_matrix, args.max_species_per_gene
        )
        gene_species_map[gene] = species_list
        gene_reasons[gene] = reason
        all_needed_species.update(species_list)

        # Check if we have any affected species
        if cand['event_type'] == 'loss':
            gene_data = loss_matrix.get(gene, {})
            affected_in_selection = [s for s in species_list if gene_data.get(s, 1) == 0]
        else:
            gene_data = dup_matrix.get(gene, {})
            affected_in_selection = [s for s in species_list if gene_data.get(s, 0) == 1]

        if not affected_in_selection:
            print(f"  WARNING: {gene} has no affected species in selection!")
            print(f"    Event type: {cand['event_type']}")
            print(f"    Selection: {species_list}")
        else:
            print(f"  {gene}: {len(affected_in_selection)} affected species selected")

    print(f"\nTotal unique species needed: {len(all_needed_species)}")

    # Get assembly info for all species
    print("\nChecking genome assemblies...")
    species_info = {}
    for species in sorted(all_needed_species):
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

    # Download genomes
    if args.skip_download:
        print("\nUsing cached genomes only...")
        genome_paths = {}
        cached_files = list(genome_dir.glob("*.fna"))
        cached_accessions = {f.stem: f for f in cached_files}
        for species, info in species_info.items():
            acc = info['accession']
            if acc in cached_accessions:
                genome_paths[species] = cached_accessions[acc]
    else:
        print("\nDownloading genomes...")
        genome_paths = {}
        for species, info in species_info.items():
            if info['ftp_path']:
                path = download_genome(info['accession'], info['ftp_path'], genome_dir)
                if path:
                    genome_paths[species] = path

    # Validate each gene - write results incrementally
    results = []

    # Determine fieldnames based on event types in candidates
    has_loss = any(c['event_type'] == 'loss' for c in candidates)
    has_dup = any(c['event_type'] == 'duplication' for c in candidates)

    fieldnames = ['gene', 'event_type', 'pattern', 'reason_flagged',
                  'n_events_total', 'primary_clade',
                  'species', 'is_affected_species', 'assembly', 'scaffold_n50',
                  'miniprot_loci', 'miniprot_raw_loci', 'miniprot_best_identity']
    if has_loss:
        fieldnames.append('toga_presence')
    if has_dup:
        fieldnames.append('toga_duplication')
    fieldnames.append('concordant')

    # Open output file for incremental writing
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        outfile.flush()

        for cand in candidates:
            gene = cand['gene']
            event_type = cand['event_type']
            reason = gene_reasons[gene]

            print(f"\n{'='*70}")
            print(f"Gene: {gene} ({event_type.upper()})")
            print(f"Reason: {reason}")
            print('='*70, flush=True)

            # Get protein sequence
            uniprot_id, protein_seq = fetch_protein_sequence(gene, session)
            if not protein_seq:
                print(f"  Could not fetch protein sequence")
                continue

            protein_fasta = genome_dir / f"{gene}.fasta"
            with open(protein_fasta, 'w') as f:
                f.write(f">{gene}\n{protein_seq}\n")
            print(f"  Protein: {uniprot_id}, {len(protein_seq)} aa", flush=True)

            # Get the appropriate TOGA matrix for this event type
            if event_type == 'loss':
                toga_data = loss_matrix
                toga_field = 'toga_presence'
            else:
                toga_data = dup_matrix
                toga_field = 'toga_duplication'

            # Validate in selected species
            species_to_check = gene_species_map.get(gene, [])

            for species in species_to_check:
                if species not in genome_paths:
                    print(f"    ✗ {species}: genome not available", flush=True)
                    continue

                genome_path = genome_paths[species]
                gff_path = genome_dir / f"{species}_{gene}.gff"

                # Run miniprot
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

                # Determine concordance
                if event_type == 'loss':
                    # Loss: presence=1 means gene present, presence=0 means lost
                    concordant = 'yes' if (mp_result['n_loci'] > 0) == (toga_value == 1) else 'no'
                    is_affected = toga_value == 0  # Lost in TOGA
                else:
                    # Duplication: 1=duplicated, 0=single copy
                    concordant = 'yes' if (mp_result['n_loci'] > 1) == (toga_value == 1) else 'no'
                    is_affected = toga_value == 1  # Duplicated in TOGA

                result = {
                    'gene': gene,
                    'event_type': event_type,
                    'pattern': cand['pattern'],
                    'reason_flagged': reason,
                    'n_events_total': cand['n_events'],
                    'primary_clade': cand['primary_clade'],
                    'species': species,
                    'is_affected_species': 'yes' if is_affected else 'no',
                    'assembly': species_info.get(species, {}).get('accession', ''),
                    'scaffold_n50': species_info.get(species, {}).get('scaffold_n50', 0),
                    'miniprot_loci': mp_result['n_loci'],
                    'miniprot_raw_loci': mp_result['n_raw_loci'],
                    'miniprot_best_identity': round(mp_result['best_identity'], 1),
                    toga_field: toga_value,
                    'concordant': concordant
                }
                results.append(result)

                # Write result immediately
                writer.writerow(result)
                outfile.flush()

                status = "✓" if concordant == 'yes' else "✗"
                affected_marker = "[AFFECTED]" if is_affected else ""
                print(f"    {status} {species} {affected_marker}: miniprot={mp_result['n_loci']} loci, TOGA={toga_value}", flush=True)

                # Clean up GFF
                if not args.keep_genomes and gff_path.exists():
                    gff_path.unlink()

    print(f"\n{'='*70}")
    print(f"Results written to {output_path}")
    print('='*70)

    # Summary
    concordant = sum(1 for r in results if r['concordant'] == 'yes')
    total = len(results)
    affected_concordant = sum(1 for r in results
                              if r['concordant'] == 'yes' and r['is_affected_species'] == 'yes')
    affected_total = sum(1 for r in results if r['is_affected_species'] == 'yes')

    print(f"\nOverall concordance: {concordant}/{total} ({100*concordant/total:.1f}%)" if total > 0 else "No results")
    print(f"Affected species concordance: {affected_concordant}/{affected_total} ({100*affected_concordant/affected_total:.1f}%)" if affected_total > 0 else "")

    # Per-gene summary
    print("\nPer-gene summary:")
    genes_seen = set()
    for r in results:
        if r['gene'] not in genes_seen:
            genes_seen.add(r['gene'])
            gene_results = [x for x in results if x['gene'] == r['gene']]
            gene_concordant = sum(1 for x in gene_results if x['concordant'] == 'yes')
            gene_affected = sum(1 for x in gene_results if x['is_affected_species'] == 'yes')
            print(f"  {r['gene']}: {gene_concordant}/{len(gene_results)} concordant, {gene_affected} affected species tested")

    return 0


if __name__ == "__main__":
    exit(main())
