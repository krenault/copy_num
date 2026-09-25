#!/usr/bin/env python3
"""Quick check of CDC42 across bat genomes."""

import subprocess
import requests
import csv
from pathlib import Path
import time

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_DELAY = 0.4

# Bat species to check
BAT_SPECIES = [
    "Myotis_brandtii",
    "Myotis_lucifugus",
    "Myotis_myotis",
    "Myotis_davidii",
    "Pteropus_giganteus",
    "Pteropus_vampyrus",
    "Eptesicus_fuscus",
]

def get_assembly_accession(species_name: str, session: requests.Session) -> str:
    """Get assembly accession for a species."""
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
        return doc.get("assemblyaccession", "")
    except:
        return None


def run_miniprot(protein_fasta: Path, genome_fasta: Path, output_gff: Path) -> bool:
    """Run miniprot."""
    try:
        cmd = ["miniprot", "-t4", "--gff", str(genome_fasta), str(protein_fasta)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=2700)
        with open(output_gff, 'w') as f:
            f.write(result.stdout)
        return True
    except Exception as e:
        print(f"    Error: {e}")
        return False


def parse_miniprot_gff(gff_path: Path) -> dict:
    """Parse miniprot GFF to count loci."""
    loci = []

    if not gff_path.exists():
        return {'n_loci': 0, 'best_identity': 0, 'loci': []}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue

            attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
            identity = float(attrs.get('Identity', 0)) * 100

            if identity >= 50:
                loci.append({
                    'chrom': parts[0],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'identity': identity
                })

    # Merge nearby loci
    loci.sort(key=lambda x: (x['chrom'], x['start']))
    merged = []
    for locus in loci:
        if not merged:
            merged.append(locus)
            continue
        last = merged[-1]
        if locus['chrom'] == last['chrom'] and locus['start'] - last['end'] < 50000:
            last['end'] = max(last['end'], locus['end'])
            last['identity'] = max(last['identity'], locus['identity'])
        else:
            merged.append(locus)

    return {
        'n_loci': len(merged),
        'best_identity': max((l['identity'] for l in merged), default=0),
        'loci': merged
    }


def main():
    genome_dir = Path("genomes")

    # Check CDC42 protein exists
    protein_fasta = genome_dir / "CDC42.fasta"
    if not protein_fasta.exists():
        print("CDC42.fasta not found, fetching...")
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": "gene:CDC42 AND organism_id:9606 AND reviewed:true", "format": "fasta", "size": 1}
        resp = requests.get(url, params=params, timeout=30)
        with open(protein_fasta, 'w') as f:
            f.write(resp.text)

    session = requests.Session()

    print("=" * 70)
    print("CDC42 duplication check across bat species")
    print("=" * 70)

    results = []

    for species in BAT_SPECIES:
        print(f"\n{species}:", flush=True)

        # Find genome
        accession = get_assembly_accession(species, session)
        if not accession:
            print(f"  No assembly found")
            continue

        genome_path = genome_dir / f"{accession}.fna"
        if not genome_path.exists():
            print(f"  Genome not cached ({accession})")
            continue

        print(f"  Assembly: {accession}", flush=True)

        # Run miniprot
        gff_path = genome_dir / f"{species}_CDC42_bat_check.gff"
        print(f"  Running miniprot...", flush=True)

        if not run_miniprot(protein_fasta, genome_path, gff_path):
            continue

        # Parse results
        result = parse_miniprot_gff(gff_path)

        print(f"  Found: {result['n_loci']} loci at {result['best_identity']:.1f}% identity")

        for i, locus in enumerate(result['loci']):
            print(f"    Locus {i+1}: {locus['chrom']}:{locus['start']}-{locus['end']} ({locus['identity']:.1f}%)")

        results.append({
            'species': species,
            'accession': accession,
            'n_loci': result['n_loci'],
            'best_identity': result['best_identity'],
            'loci': result['loci']
        })

        # Cleanup
        if gff_path.exists():
            gff_path.unlink()

    # Summary
    print("\n" + "=" * 70)
    print("Summary: CDC42 copy number in bats")
    print("=" * 70)

    for r in results:
        dup_status = "DUPLICATED" if r['n_loci'] > 1 else "single copy"
        print(f"  {r['species']}: {r['n_loci']} copies ({dup_status})")

    # Save to CSV
    output_path = Path("results/gene_validations/CDC42_bat_duplication.csv")
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['species', 'accession', 'n_loci', 'best_identity', 'loci_details'])
        for r in results:
            loci_str = "; ".join([f"{l['chrom']}:{l['start']}-{l['end']}" for l in r['loci']])
            writer.writerow([r['species'], r['accession'], r['n_loci'], r['best_identity'], loci_str])

    print(f"\nResults saved to: {output_path}")


if __name__ == "__main__":
    main()
