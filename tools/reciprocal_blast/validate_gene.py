#!/usr/bin/env python3
"""
Reciprocal BLAST (API-only) pipeline using UniProt/Ensembl human references.

Workflow:
1) Fetch human reference sequence from UniProt (protein) or Ensembl (CDS).
2) Forward search (BLAST API) against target taxon (tblastn for protein, blastn for CDS).
3) Reciprocal search of top hit sequence back to human (txid9606).
4) Accept homolog if the original gene is the top reciprocal hit.

Species without NCBI sequence data are marked no_genome and excluded from loss counts.
"""

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import List, Optional, Tuple
from xml.etree import ElementTree

import requests


SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output" / "copy_num_BLAT_check"

# Default tree file (RAxML mammalian phylogeny)
DEFAULT_TREE_FILE = DATA_DIR / "RAxML_bipartitions.result_FIN4_raw_rooted_wBoots_4098mam1out_OK.newick"

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
ENSEMBL_URL = "https://rest.ensembl.org"
UNIPROT_URL = "https://rest.uniprot.org"
NCBI_TAXONOMY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Rate limiting for NCBI
NCBI_DELAY = 0.4  # seconds between requests (max 3/sec without API key)


def ncbi_request(session: requests.Session, url: str, params: dict,
                 retries: int = 3, timeout: int = 30) -> Optional[dict]:
    """Make NCBI API request with retry logic for rate limiting."""
    for attempt in range(retries):
        try:
            time.sleep(NCBI_DELAY)  # Rate limit
            resp = session.get(url, params=params, timeout=timeout)

            if resp.status_code == 429:
                wait_time = (attempt + 1) * 10  # Exponential backoff
                print(f"    Rate limited (429), waiting {wait_time}s...")
                time.sleep(wait_time)
                continue

            resp.raise_for_status()
            return resp.json()

        except requests.exceptions.HTTPError as e:
            if "429" in str(e) and attempt < retries - 1:
                wait_time = (attempt + 1) * 10
                print(f"    Rate limited, waiting {wait_time}s...")
                time.sleep(wait_time)
                continue
            raise
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            if attempt < retries - 1:
                wait_time = (attempt + 1) * 5
                print(f"    Connection error, retrying in {wait_time}s...")
                time.sleep(wait_time)
                continue
            raise

    return None


EXCLUDE_PATTERNS = [
    r"^OR\d+",
    r"^ORF\d+",
    r"^ZNF",
    r"^HIST",
    r"^HLA",
]

# Major mammalian orders for outgroup sampling
MAMMAL_ORDERS = [
    "PRIMATES",
    "RODENTIA",
    "CARNIVORA",
    "ARTIODACTYLA",
    "CETACEA",
    "CHIROPTERA",
    "PERISSODACTYLA",
    "LAGOMORPHA",
    "EULIPOTYPHLA",
    "AFROTHERIA",  # Includes Proboscidea, Sirenia, etc.
    "XENARTHRA",
    "MARSUPIALIA",
]


def extract_order_from_species(species_name: str) -> Optional[str]:
    """Extract taxonomic order from species name (e.g., 'Homo_sapiens_HOMINIDAE_PRIMATES' -> 'PRIMATES')."""
    parts = species_name.split("_")
    if len(parts) >= 4:
        # Last part is typically the order
        potential_order = parts[-1].upper()
        # Check if it looks like an order (all caps, no numbers)
        if potential_order.isalpha() and potential_order.isupper():
            return potential_order
    return None


def get_genome_quality_metrics(taxid: int, session: requests.Session) -> dict:
    """Fetch genome assembly quality metrics from NCBI Assembly database."""
    metrics = {
        "assembly_accession": "",
        "assembly_level": "",
        "contig_count": "",
        "contig_n50": "",
        "scaffold_count": "",
        "scaffold_n50": "",
        "total_length": "",
    }

    try:
        # Search for assembly
        url = f"{NCBI_TAXONOMY_URL}/esearch.fcgi"
        params = {
            "db": "assembly",
            "term": f"txid{taxid}[ORGN] AND (latest[filter] OR \"representative genome\"[filter])",
            "retmode": "json",
            "retmax": 1,
        }

        data = ncbi_request(session, url, params)
        if not data:
            return metrics

        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            return metrics

        assembly_id = ids[0]

        # Fetch assembly summary
        summary_url = f"{NCBI_TAXONOMY_URL}/esummary.fcgi"
        summary_params = {
            "db": "assembly",
            "id": assembly_id,
            "retmode": "json",
        }

        summary = ncbi_request(session, summary_url, summary_params)
        if not summary:
            return metrics

        doc = summary.get("result", {}).get(assembly_id, {})
        if doc:
            metrics["assembly_accession"] = doc.get("assemblyaccession", "")
            metrics["assembly_level"] = doc.get("assemblystatus", "")
            metrics["contig_count"] = doc.get("contign50", "")  # This is actually N50
            metrics["contig_n50"] = doc.get("contign50", "")
            metrics["scaffold_n50"] = doc.get("scaffoldn50", "")
            metrics["total_length"] = doc.get("totallength", "")

            # Try to get contig/scaffold counts from stats
            stats = doc.get("assemblystatistics", {})
            if isinstance(stats, dict):
                metrics["contig_count"] = stats.get("contigcount", "")
                metrics["scaffold_count"] = stats.get("scaffoldcount", "")

    except Exception as e:
        print(f"    Warning: Could not fetch assembly metrics: {e}")

    return metrics


def load_neighbors_from_file(file_path: Path, n_neighbors: int, focal_species: str) -> List[str]:
    neighbors = []
    with open(file_path) as f:
        for line in f:
            species = line.strip()
            if not species:
                continue
            if species.replace(" ", "_").lower() == focal_species.replace(" ", "_").lower():
                continue
            neighbors.append(species)
            if len(neighbors) >= n_neighbors:
                break
    return neighbors


def load_neighbors_from_tree(tree_file: Path, focal_species: str, n_neighbors: int) -> List[str]:
    """Load N closest neighbors from tree (does not check NCBI availability)."""
    try:
        from Bio import Phylo  # type: ignore
    except Exception as exc:
        raise SystemExit("BioPython is required for --tree-file neighbor selection.") from exc

    tree = Phylo.read(str(tree_file), "newick")
    terminals = [clade for clade in tree.get_terminals() if clade.name]

    focal_key = focal_species.replace(" ", "_")
    focal_clade = None
    for clade in terminals:
        if clade.name and clade.name.lower().startswith(focal_key.lower()):
            focal_clade = clade
            break

    if not focal_clade:
        raise SystemExit(f"Focal species '{focal_species}' not found in tree")

    distances = []
    for clade in terminals:
        if clade.name and clade.name != focal_clade.name:
            try:
                dist = tree.distance(focal_clade, clade)
                distances.append((clade.name, dist))
            except Exception:
                continue

    distances.sort(key=lambda x: x[1])
    return [name for name, _ in distances[:n_neighbors]]


def get_all_neighbors_sorted(tree_file: Path, focal_species: str) -> List[str]:
    """Get ALL neighbors from tree, sorted by phylogenetic distance."""
    try:
        from Bio import Phylo  # type: ignore
    except Exception as exc:
        raise SystemExit("BioPython is required for tree neighbor selection.") from exc

    tree = Phylo.read(str(tree_file), "newick")
    terminals = [clade for clade in tree.get_terminals() if clade.name]

    focal_key = focal_species.replace(" ", "_")
    focal_clade = None
    for clade in terminals:
        if clade.name and clade.name.lower().startswith(focal_key.lower()):
            focal_clade = clade
            break

    if not focal_clade:
        raise SystemExit(f"Focal species '{focal_species}' not found in tree")

    distances = []
    for clade in terminals:
        if clade.name and clade.name != focal_clade.name:
            try:
                dist = tree.distance(focal_clade, clade)
                distances.append((clade.name, dist))
            except Exception:
                continue

    distances.sort(key=lambda x: x[1])
    return [name for name, _ in distances]


def find_neighbors_with_ncbi_data(
    tree_file: Path,
    focal_species: str,
    n_neighbors: int,
    session: requests.Session,
    include_focal: bool = True,
    max_search: int = 100,
    min_sequences: int = 10000,
) -> Tuple[List[str], List[str]]:
    """
    Find N neighbors that have NCBI genome data (sufficient sequences).

    Returns:
        Tuple of (species_with_data, species_without_data)
    """
    print(f"\nSearching for {n_neighbors} species with NCBI genome (>={min_sequences} sequences)...")

    # Get all neighbors sorted by distance
    all_neighbors = get_all_neighbors_sorted(tree_file, focal_species)

    species_with_data = []
    species_without_data = []
    checked = 0

    # Optionally check focal species first
    if include_focal:
        focal_taxid = get_ncbi_taxid(focal_species, session)
        if focal_taxid:
            count = ncbi_nuccore_count(focal_taxid, session)
            if count >= min_sequences:
                print(f"  + {focal_species}: {count:,} sequences (GENOME)")
                species_with_data.append(focal_species)
            elif count > 0:
                print(f"  - {focal_species}: {count:,} sequences (no genome)")
                species_without_data.append(focal_species)
            else:
                print(f"  - {focal_species}: no NCBI data")
                species_without_data.append(focal_species)
        else:
            print(f"  - {focal_species}: no taxid found")
            species_without_data.append(focal_species)

    # Search through neighbors until we find N with data
    for neighbor in all_neighbors:
        if len(species_with_data) >= n_neighbors + (1 if include_focal else 0):
            break

        if checked >= max_search:
            print(f"  Reached max search limit ({max_search})")
            break

        checked += 1

        # Clean species name for lookup
        clean_name = neighbor.split("_")[0] + "_" + neighbor.split("_")[1] if "_" in neighbor else neighbor

        taxid = get_ncbi_taxid(clean_name, session)
        if not taxid:
            print(f"  - {clean_name}: no taxid")
            species_without_data.append(neighbor)
            continue

        count = ncbi_nuccore_count(taxid, session)
        if count >= min_sequences:
            print(f"  + {clean_name}: {count:,} sequences (GENOME)")
            species_with_data.append(neighbor)
        elif count > 0:
            print(f"  - {clean_name}: {count:,} sequences (no genome)")
            species_without_data.append(neighbor)
        else:
            print(f"  - {clean_name}: no NCBI data")
            species_without_data.append(neighbor)

    found = len(species_with_data) - (1 if include_focal and species_with_data else 0)
    print(f"\nFound {found} neighbors with genome (checked {checked} species)")

    return species_with_data, species_without_data


def find_outgroups_with_ncbi_data(
    tree_file: Path,
    focal_species: str,
    n_outgroups: int,
    session: requests.Session,
    exclude_species: List[str] = None,
    max_search: int = 100,
    min_sequences: int = 10000,
) -> Tuple[List[str], List[str]]:
    """
    Find N outgroup species (distant from focal) that have NCBI genome data.
    Samples from the most distant species in the tree.

    Returns:
        Tuple of (species_with_data, species_without_data)
    """
    print(f"\nSearching for {n_outgroups} outgroup species with genome (>={min_sequences} sequences)...")

    # Get all neighbors sorted by distance (closest first)
    all_neighbors = get_all_neighbors_sorted(tree_file, focal_species)

    # Reverse to get most distant first
    all_neighbors = list(reversed(all_neighbors))

    # Exclude species already in focal clade sample
    if exclude_species:
        exclude_set = set(s.replace(" ", "_").lower() for s in exclude_species)
        all_neighbors = [n for n in all_neighbors
                        if n.replace(" ", "_").lower() not in exclude_set]

    species_with_data = []
    species_without_data = []
    checked = 0

    for neighbor in all_neighbors:
        if len(species_with_data) >= n_outgroups:
            break

        if checked >= max_search:
            print(f"  Reached max search limit ({max_search})")
            break

        checked += 1

        # Clean species name for lookup
        clean_name = neighbor.split("_")[0] + "_" + neighbor.split("_")[1] if "_" in neighbor else neighbor

        taxid = get_ncbi_taxid(clean_name, session)
        if not taxid:
            print(f"  - {clean_name}: no taxid (outgroup)")
            species_without_data.append(neighbor)
            continue

        count = ncbi_nuccore_count(taxid, session)
        if count >= min_sequences:
            print(f"  + {clean_name}: {count:,} sequences (GENOME, outgroup)")
            species_with_data.append(neighbor)
        elif count > 0:
            print(f"  - {clean_name}: {count:,} sequences (no genome, outgroup)")
            species_without_data.append(neighbor)
        else:
            print(f"  - {clean_name}: no NCBI data (outgroup)")
            species_without_data.append(neighbor)

    print(f"\nFound {len(species_with_data)} outgroups with genome (checked {checked} species)")

    return species_with_data, species_without_data


def find_outgroups_by_order(
    tree_file: Path,
    focal_species: str,
    n_per_order: int,
    orders_to_sample: List[str],
    session: requests.Session,
    exclude_species: List[str] = None,
    max_search_per_order: int = 50,
    min_sequences: int = 10000,
) -> Tuple[List[str], dict]:
    """
    Find N species from each specified mammalian order that have NCBI genome data.

    Returns:
        Tuple of (species_with_data, order_counts_dict)
    """
    try:
        from Bio import Phylo
    except ImportError:
        raise SystemExit("BioPython is required for tree-based sampling.")

    tree = Phylo.read(str(tree_file), "newick")
    terminals = [clade.name for clade in tree.get_terminals() if clade.name]

    # Get focal species order to exclude it
    focal_key = focal_species.replace(" ", "_")
    focal_order = None
    for term in terminals:
        if term.lower().startswith(focal_key.lower()):
            focal_order = extract_order_from_species(term)
            break

    print(f"\n=== Sampling outgroups by order ===")
    print(f"Focal species: {focal_species} (Order: {focal_order})")
    print(f"Sampling {n_per_order} species per order with genome (>={min_sequences:,} sequences)")
    print(f"Orders to sample: {', '.join(orders_to_sample)}")

    # Exclude focal clade species
    exclude_set = set()
    if exclude_species:
        exclude_set = set(s.replace(" ", "_").lower() for s in exclude_species)

    # Group tree species by order
    species_by_order = {}
    for term in terminals:
        order = extract_order_from_species(term)
        if order and order != focal_order:  # Exclude focal order
            if order not in species_by_order:
                species_by_order[order] = []
            if term.replace(" ", "_").lower() not in exclude_set:
                species_by_order[order].append(term)

    all_species_with_data = []
    order_counts = {}

    for order in orders_to_sample:
        order_upper = order.upper()
        if order_upper not in species_by_order:
            print(f"\n  {order_upper}: not found in tree (or same as focal)")
            order_counts[order_upper] = 0
            continue

        candidates = species_by_order[order_upper]
        print(f"\n  {order_upper}: {len(candidates)} candidates in tree")

        found_for_order = []
        checked = 0

        for species in candidates:
            if len(found_for_order) >= n_per_order:
                break
            if checked >= max_search_per_order:
                print(f"    Reached search limit ({max_search_per_order})")
                break

            checked += 1

            # Clean species name
            clean_name = species.split("_")[0] + "_" + species.split("_")[1] if "_" in species else species

            taxid = get_ncbi_taxid(clean_name, session)
            if not taxid:
                continue

            count = ncbi_nuccore_count(taxid, session)
            if count >= min_sequences:
                print(f"    + {clean_name}: {count:,} sequences (GENOME)")
                found_for_order.append(species)
                all_species_with_data.append(species)
            elif count > 1000:  # Show near-misses
                print(f"    - {clean_name}: {count:,} sequences (below threshold)")

        order_counts[order_upper] = len(found_for_order)
        print(f"    Found {len(found_for_order)}/{n_per_order} for {order_upper}")

    print(f"\n=== Total outgroups: {len(all_species_with_data)} species across {len([o for o in order_counts if order_counts[o] > 0])} orders ===")

    return all_species_with_data, order_counts


def get_ncbi_taxid(species_name: str, session: requests.Session) -> Optional[int]:
    """Look up NCBI taxonomy ID for a species."""
    clean_name = species_name.replace("_", " ")
    parts = clean_name.split()
    if len(parts) > 2:
        clean_name = " ".join(parts[:2])

    url = f"{NCBI_TAXONOMY_URL}/esearch.fcgi"
    params = {
        "db": "taxonomy",
        "term": clean_name,
        "retmode": "json",
    }

    try:
        data = ncbi_request(session, url, params)
        if data:
            ids = data.get("esearchresult", {}).get("idlist", [])
            if ids:
                return int(ids[0])
    except Exception as e:
        print(f"  Warning: Could not get taxid for {species_name}: {e}")

    return None


def ncbi_nuccore_count(taxid: int, session: requests.Session) -> int:
    """Return count of nuccore entries for a taxid."""
    url = f"{NCBI_TAXONOMY_URL}/esearch.fcgi"
    params = {
        "db": "nuccore",
        "term": f"txid{taxid}[ORGN]",
        "retmode": "json",
        "retmax": 0,
    }

    try:
        data = ncbi_request(session, url, params)
        if data:
            count = int(data.get("esearchresult", {}).get("count", 0))
            return count
    except Exception as e:
        print(f"  Warning: Could not check nuccore count for txid{taxid}: {e}")

    return 0


def get_focal_sequence(gene_symbol: str, focal_taxid: int, session: requests.Session) -> dict:
    """Fetch a focal species nucleotide sequence for a gene symbol from NCBI nuccore."""
    url = f"{NCBI_TAXONOMY_URL}/esearch.fcgi"
    params = {
        "db": "nuccore",
        "term": f"{gene_symbol}[Gene] AND txid{focal_taxid}[ORGN]",
        "retmode": "json",
        "retmax": 1,
    }

    try:
        resp = session.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            accession = ids[0]
            seq = fetch_ncbi_sequence(accession, session, db="nuccore")
            if seq:
                print(f"  Found focal sequence accession: {accession} ({len(seq)} bp)")
                print(f"  Debug: focal sequence prefix: {seq[:30]}")
                return {
                    "id": accession,
                    "seq": seq,
                    "source": "focal_ncbi",
                }
    except Exception as e:
        print(f"  Warning: Could not fetch focal sequence for {gene_symbol}: {e}")

    return {"id": "", "seq": "", "source": "focal_ncbi"}


def get_gene_id_from_symbol(gene_symbol: str, session: requests.Session) -> Optional[str]:
    """Look up Ensembl gene ID from gene symbol (human)."""
    url = f"{ENSEMBL_URL}/xrefs/symbol/homo_sapiens/{gene_symbol}"
    headers = {"Content-Type": "application/json"}

    try:
        resp = session.get(url, headers=headers, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        for entry in data:
            if entry.get("type") == "gene":
                return entry.get("id")

        if data:
            return data[0].get("id")

    except Exception as e:
        print(f"Error looking up gene {gene_symbol}: {e}")

    return None


def get_uniprot_sequence(gene_symbol: str, session: requests.Session) -> dict:
    """Get canonical protein sequence from UniProt for human gene."""
    # Search UniProt for human gene
    url = f"{UNIPROT_URL}/uniprotkb/search"
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 1,
    }

    try:
        resp = session.get(url, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        results = data.get("results", [])
        if results:
            entry = results[0]
            seq = entry.get("sequence", {}).get("value", "")
            accession = entry.get("primaryAccession", "")
            if seq:
                print(f"  Debug: UniProt sequence prefix: {seq[:30]}")
            return {
                "id": accession,
                "seq": seq,
                "source": "uniprot",
                "gene_symbol": gene_symbol,
            }
    except Exception as e:
        print(f"  UniProt lookup failed: {e}")

    return {"id": "", "seq": "", "source": "uniprot"}


def get_ensembl_sequence(gene_id: str, session: requests.Session, seq_type: str = "cds") -> dict:
    """Get gene sequence from Ensembl (human)."""
    url = f"{ENSEMBL_URL}/sequence/id/{gene_id}"
    params = {"type": seq_type, "multiple_sequences": "1"}
    headers = {"Content-Type": "application/json"}

    try:
        resp = session.get(url, headers=headers, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()

        if isinstance(data, list):
            if data:
                longest = max(data, key=lambda x: len(x.get("seq", "")))
                return {
                    "id": longest.get("id", gene_id),
                    "seq": longest.get("seq", ""),
                    "source": "ensembl",
                }
            return {"id": gene_id, "seq": "", "source": "ensembl"}
        else:
            return {
                "id": data.get("id", gene_id),
                "seq": data.get("seq", ""),
                "source": "ensembl",
            }
    except Exception as e:
        print(f"Error fetching sequence for {gene_id}: {e}")
        return {"id": gene_id, "seq": "", "source": "ensembl"}


def get_human_sequence(gene_symbol: str, session: requests.Session,
                       use_protein: bool = False,
                       reference_source: str = "auto") -> dict:
    """
    Get human reference sequence for a gene.

    Tries UniProt first (for protein), then Ensembl (for CDS/cDNA).
    Human is used as reference since TOGA uses human as the reference genome.
    """
    print(f"\n[Fetching human reference sequence for {gene_symbol}]")

    reference_source = reference_source.lower()
    if reference_source not in {"auto", "uniprot", "ensembl"}:
        reference_source = "auto"

    if use_protein or reference_source == "uniprot":
        # Try UniProt for protein sequence
        print("  Trying UniProt (protein)...")
        result = get_uniprot_sequence(gene_symbol, session)
        if result["seq"]:
            print(f"  Found: UniProt {result['id']} ({len(result['seq'])} aa)")
            return result

        if reference_source == "uniprot":
            print("  UniProt requested but not found; falling back to Ensembl.")

    # Try Ensembl for CDS
    print("  Trying Ensembl (CDS)...")
    gene_id = get_gene_id_from_symbol(gene_symbol, session)
    if gene_id:
        result = get_ensembl_sequence(gene_id, session, seq_type="cds")
        if result["seq"]:
            result["gene_id"] = gene_id
            print(f"  Found: Ensembl {gene_id} ({len(result['seq'])} bp)")
            return result

    print(f"  ERROR: Could not find sequence for {gene_symbol}")
    return {"id": "", "seq": "", "source": "none"}


def fetch_ncbi_sequence(accession: str, session: requests.Session, db: str = "nuccore") -> str:
    """Fetch a nucleotide/protein sequence from NCBI by accession."""
    url = f"{NCBI_TAXONOMY_URL}/efetch.fcgi"
    params = {
        "db": db,
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }

    try:
        resp = session.get(url, params=params, timeout=60)
        resp.raise_for_status()
        fasta = resp.text
        if not fasta:
            return ""
        lines = [line.strip() for line in fasta.splitlines() if line.strip()]
        if len(lines) <= 1:
            return ""
        return "".join(lines[1:])
    except Exception as e:
        print(f"  Warning: could not fetch sequence for {accession}: {e}")
        return ""


def reciprocal_blast(top_hit_accession: str,
                     gene_symbol: str,
                     session: requests.Session,
                     target_taxid: int,
                     max_seq_len: int = 8000) -> Tuple[str, dict]:
    """Run reciprocal BLAST of top hit back to a target taxon; returns status and details."""
    if not top_hit_accession:
        return "no_accession", {}

    sequence = fetch_ncbi_sequence(top_hit_accession, session, db="nuccore")
    if not sequence:
        return "no_sequence", {}

    if len(sequence) > max_seq_len:
        return "skipped_too_long", {"sequence_length": len(sequence)}

    entrez_query = f"txid{target_taxid}[ORGN]"
    rid = submit_blast(sequence, "nt", "blastn", entrez_query, session)
    if not rid:
        return "blast_error", {}

    hits = wait_and_get_blast(rid, session)
    if not hits:
        return "no_hits", {}

    top_hit = hits[0]
    gene_lower = gene_symbol.lower()
    matching = [h for h in hits if gene_lower in h.get("definition", "").lower()]
    status = "confirmed" if matching else "no_match"
    return status, {
        "top_hit_accession": top_hit.get("accession", ""),
        "top_hit_definition": top_hit.get("definition", ""),
        "top_hit_identity": top_hit.get("pct_identity", 0),
        "top_hit_evalue": top_hit.get("evalue", ""),
    }


def submit_blast(sequence: str, database: str, program: str,
                 entrez_query: str, session: requests.Session,
                 expect: float = 0.001) -> Optional[str]:
    """Submit BLAST job and return RID."""
    print(f"  Debug: submitting BLAST ({program} vs {database}), query_len={len(sequence)}, entrez='{entrez_query}'")
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "EXPECT": expect,
        "FORMAT_TYPE": "XML",
        "HITLIST_SIZE": 50,
    }

    if entrez_query:
        params["ENTREZ_QUERY"] = entrez_query

    try:
        resp = session.post(BLAST_URL, data=params, timeout=60)
        resp.raise_for_status()

        match = re.search(r"RID = (\w+)", resp.text)
        if match:
            return match.group(1)
        print("  Debug: BLAST response did not include RID")
    except Exception as e:
        print(f"  BLAST submission error: {e}")
        try:
            snippet = resp.text[:300] if resp is not None else ""
            if snippet:
                print(f"  Debug: BLAST response snippet: {snippet}")
        except Exception:
            pass

    return None


def check_blast_status(rid: str, session: requests.Session, retries: int = 3) -> str:
    """Check BLAST job status with retry logic."""
    params = {
        "CMD": "Get",
        "FORMAT_OBJECT": "SearchInfo",
        "RID": rid,
    }

    for attempt in range(retries):
        try:
            resp = session.get(BLAST_URL, params=params, timeout=60)
            resp.raise_for_status()

            if "Status=WAITING" in resp.text:
                return "WAITING"
            elif "Status=READY" in resp.text:
                return "READY"
            elif "Status=FAILED" in resp.text:
                return "FAILED"
            return "UNKNOWN"
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            if attempt < retries - 1:
                wait_time = (attempt + 1) * 5
                print(f"  Connection error, retrying in {wait_time}s... ({e.__class__.__name__})")
                time.sleep(wait_time)
            else:
                print(f"  Connection failed after {retries} attempts")
                return "ERROR"
        except Exception as e:
            print(f"  Unexpected error checking BLAST status: {e}")
            return "ERROR"

    return "UNKNOWN"


def get_blast_results(rid: str, session: requests.Session) -> str:
    """Get BLAST results XML."""
    params = {
        "CMD": "Get",
        "FORMAT_TYPE": "XML",
        "RID": rid,
    }

    resp = session.get(BLAST_URL, params=params, timeout=120)
    resp.raise_for_status()
    return resp.text


def parse_blast_xml(xml_text: str) -> list:
    """Parse BLAST XML and return hits."""
    hits = []

    try:
        root = ElementTree.fromstring(xml_text)

        for hit in root.findall(".//Hit"):
            hit_def = hit.find("Hit_def").text if hit.find("Hit_def") is not None else ""
            hit_accession = hit.find("Hit_accession").text if hit.find("Hit_accession") is not None else ""
            hit_len = int(hit.find("Hit_len").text) if hit.find("Hit_len") is not None else 0

            best_hsp = None
            best_score = 0

            for hsp in hit.findall(".//Hsp"):
                score = float(hsp.find("Hsp_bit-score").text) if hsp.find("Hsp_bit-score") is not None else 0
                if score > best_score:
                    best_score = score
                    best_hsp = {
                        "bit_score": score,
                        "evalue": float(hsp.find("Hsp_evalue").text) if hsp.find("Hsp_evalue") is not None else 999,
                        "identity": int(hsp.find("Hsp_identity").text) if hsp.find("Hsp_identity") is not None else 0,
                        "align_len": int(hsp.find("Hsp_align-len").text) if hsp.find("Hsp_align-len") is not None else 0,
                    }

            if best_hsp:
                pct_identity = (best_hsp["identity"] / best_hsp["align_len"] * 100) if best_hsp["align_len"] > 0 else 0
                hits.append({
                    "accession": hit_accession,
                    "definition": hit_def,
                    "length": hit_len,
                    "bit_score": best_hsp["bit_score"],
                    "evalue": best_hsp["evalue"],
                    "pct_identity": pct_identity,
                    "align_len": best_hsp["align_len"],
                })

    except ElementTree.ParseError as e:
        print(f"  XML parse error: {e}")

    return sorted(hits, key=lambda x: x["bit_score"], reverse=True)


def wait_and_get_blast(rid: str, session: requests.Session,
                       max_wait: int = 300, poll_interval: int = 15) -> list:
    """Wait for BLAST to complete and return parsed results."""
    elapsed = 0
    while elapsed < max_wait:
        status = check_blast_status(rid, session)
        if status == "READY":
            xml_results = get_blast_results(rid, session)
            return parse_blast_xml(xml_results)
        elif status == "FAILED":
            print(f"  BLAST failed")
            return []
        elif status == "ERROR":
            print(f"  BLAST connection error, skipping...")
            return []

        print(f"  Waiting... ({elapsed}s)")
        time.sleep(poll_interval)
        elapsed += poll_interval

    print(f"  BLAST timed out")
    return []


def blast_species(sequence: str, species_name: str, taxid: int,
                  session: requests.Session, gene_symbol: str,
                  use_protein: bool = False,
                  reciprocal: bool = True,
                  reciprocal_max_len: int = 8000) -> dict:
    """BLAST reference sequence against a species and return results."""
    clean_name = species_name.replace("_", " ").split()[0:2]
    clean_name = " ".join(clean_name)

    print(f"\n--- {clean_name} (taxid:{taxid}) ---")

    entrez_query = f"txid{taxid}[ORGN]"

    # Choose program based on sequence type
    if use_protein:
        program = "tblastn"  # Protein query vs translated nucleotide DB
    else:
        program = "blastn"

    query_type = "protein" if use_protein else "nucleotide"

    # Try multiple databases: nt first, then refseq_genomes, then wgs
    # WGS and refseq_genomes need organism name filter instead of txid
    databases_to_try = [
        ("nt", entrez_query),                                      # nt with txid filter
        ("refseq_genomes", f'"{clean_name}"[ORGN]'),              # refseq with species name
        ("wgs", f'"{clean_name}"[ORGN]'),                          # wgs with species name
        ("wgs", ""),                                               # wgs with no filter (last resort)
    ]
    hits = None

    for database, db_entrez in databases_to_try:
        print(f"  Debug: query_type={query_type}, program={program}, db={database}, seq_len={len(sequence)}")
        if db_entrez:
            print(f"  Using entrez filter: {db_entrez}")

        rid = submit_blast(sequence, database, program, db_entrez, session)
        if not rid:
            print(f"  Could not submit BLAST to {database}, trying next...")
            continue

        print(f"  RID: {rid}")
        hits = wait_and_get_blast(rid, session)

        if hits:
            print(f"  Found {len(hits)} hits in {database}")
            break
        else:
            print(f"  No hits in {database}, trying next database...")

    if hits:
        top = hits[0]
        total_aligned = sum(h.get("align_len", 0) for h in hits)
        query_len = len(sequence) if sequence else 1
        ecnc = total_aligned / query_len
        n_hits = len(hits)
        duplicated_call = n_hits >= 2 and ecnc >= 1.5
        result = {
            "status": "found",
            "n_hits": n_hits,
            "ecnc": round(ecnc, 3),
            "duplicated_call": duplicated_call,
            "top_hit_accession": top["accession"],
            "top_hit_definition": top["definition"],
            "top_hit_identity": top["pct_identity"],
            "top_hit_evalue": top["evalue"],
        }
        if reciprocal:
            status, details = reciprocal_blast(
                result["top_hit_accession"],
                gene_symbol,
                session,
                target_taxid=9606,
                max_seq_len=reciprocal_max_len,
            )
            result["reciprocal_status"] = status
            result["reciprocal_details"] = details
        return result

    if ncbi_nuccore_count(taxid, session) == 0:
        return {"status": "no_genome"}

    return {"status": "no_hits"}


def validate_gene_single(gene_symbol: str, target_species: str, reference_source: str,
                         use_protein: bool, reciprocal_max_len: int = 8000,
                         fetch_genome_metrics: bool = True) -> dict:
    """Run reciprocal BLAST for one gene and one target species."""
    session = requests.Session()

    taxid = get_ncbi_taxid(target_species, session)
    if not taxid:
        return {
            "gene_symbol": gene_symbol,
            "target_species": target_species,
            "status": "no_taxid",
        }

    # Fetch genome quality metrics
    genome_metrics = {}
    if fetch_genome_metrics and taxid:
        genome_metrics = get_genome_quality_metrics(taxid, session)
        
    nuccore_count = ncbi_nuccore_count(taxid, session)
    if nuccore_count == 0:
        result = {
            "gene_symbol": gene_symbol,
            "target_species": target_species,
            "status": "no_genome",
            "nuccore_count": nuccore_count,
        }
        result.update(genome_metrics)
        return result

    seq_data = get_human_sequence(
        gene_symbol,
        session,
        use_protein=use_protein,
        reference_source=reference_source,
    )

    if not seq_data["seq"]:
        return {
            "gene_symbol": gene_symbol,
            "target_species": target_species,
            "status": "no_reference",
        }

    result = blast_species(
        seq_data["seq"],
        target_species,
        taxid,
        session,
        gene_symbol,
        use_protein=use_protein,
        reciprocal=True,
        reciprocal_max_len=reciprocal_max_len,
    )

    result.update({
        "gene_symbol": gene_symbol,
        "target_species": target_species,
        "reference_source": seq_data["source"],
        "reference_id": seq_data["id"],
        "reference_length": len(seq_data["seq"]),
        "nuccore_count": nuccore_count,
    })
    result.update(genome_metrics)

    return result


def validate_gene_multi(gene_symbol: str, species_list: List[str], reference_source: str,
                        use_protein: bool, reciprocal_max_len: int = 8000) -> List[dict]:
    results = []
    for i, species in enumerate(species_list):
        print(f"\n[{i+1}/{len(species_list)}] Processing {species}...")
        try:
            result = validate_gene_single(
                gene_symbol=gene_symbol,
                target_species=species,
                reference_source=reference_source,
                use_protein=use_protein,
                reciprocal_max_len=reciprocal_max_len,
            )
            results.append(result)
        except Exception as e:
            print(f"  ERROR processing {species}: {e}")
            results.append({
                "gene_symbol": gene_symbol,
                "target_species": species,
                "status": "error",
                "error_message": str(e),
            })
        time.sleep(1)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Reciprocal BLAST pipeline (API-only)",
    )
    parser.add_argument("--gene", required=True,
                        help="Human gene symbol (e.g., NKG7)")
    parser.add_argument("--target-species", default=None,
                        help="Target species name (e.g., Myotis davidii)")
    parser.add_argument("--focal-species", default=None,
                        help="Focal species for neighbor selection")
    parser.add_argument("--n-neighbors", type=int, default=0,
                        help="Number of neighbors to include (uses RAxML tree by default)")
    parser.add_argument("--require-ncbi", action="store_true",
                        help="Only include species with NCBI data (search until N found)")
    parser.add_argument("--outgroup-samples", type=int, default=0,
                        help="Number of distant outgroup species to sample (verifies clade-specificity)")
    parser.add_argument("--outgroups-per-order", type=int, default=0,
                        help="Sample N species from each major mammalian order (more thorough than --outgroup-samples)")
    parser.add_argument("--orders", type=str, default=None,
                        help="Comma-separated list of orders to sample (default: PRIMATES,RODENTIA,CARNIVORA,ARTIODACTYLA,CETACEA)")
    parser.add_argument("--min-sequences", type=int, default=10000,
                        help="Minimum NCBI sequences to consider species has genome (default: 10000)")
    parser.add_argument("--max-search", type=int, default=100,
                        help="Max species to check when using --require-ncbi (default: 100)")
    parser.add_argument("--neighbors-file", type=Path, default=None,
                        help="File with one neighbor species per line (optional)")
    parser.add_argument("--tree-file", type=Path, default=None,
                        help=f"Newick tree file (default: {DEFAULT_TREE_FILE.name})")
    parser.add_argument("--reference-source", choices=["auto", "uniprot", "ensembl"],
                        default="auto",
                        help="Human reference source: auto/uniprot/ensembl (default: auto)")
    parser.add_argument("--protein", action="store_true",
                        help="Use protein query (tblastn). Default uses CDS/blastn.")
    parser.add_argument("--reciprocal-max-length", type=int, default=8000,
                        help="Max length for reciprocal sequence (default: 8000)")
    parser.add_argument("--output", "-o", type=Path, default=None,
                        help="Output CSV file (default: auto-generated)")
    parser.add_argument("--write-json", action="store_true",
                        help="Also write full JSON output for debugging")

    args = parser.parse_args()

    if args.reference_source == "uniprot" and not args.protein:
        print("Note: UniProt reference is protein; enabling --protein mode.")
        args.protein = True

    gene_symbol = args.gene
    for pat in EXCLUDE_PATTERNS:
        if re.match(pat, gene_symbol, re.IGNORECASE):
            raise SystemExit(f"Gene {gene_symbol} excluded by pattern {pat}.")

    if args.target_species:
        species_list = [args.target_species]
    else:
        if not args.focal_species:
            raise SystemExit("Provide --target-species or --focal-species.")

        # Use provided tree file or default RAxML tree
        tree_file = args.tree_file if args.tree_file else DEFAULT_TREE_FILE
        if not tree_file.exists():
            raise SystemExit(f"Tree file not found: {tree_file}")
        print(f"Using tree: {tree_file.name}")

        if args.n_neighbors > 0:
            if args.neighbors_file:
                neighbors = load_neighbors_from_file(args.neighbors_file, args.n_neighbors, args.focal_species)
                species_list = [args.focal_species] + neighbors
            elif args.require_ncbi:
                # Find N neighbors that have NCBI genome data
                session = requests.Session()
                species_with_data, _ = find_neighbors_with_ncbi_data(
                    tree_file=tree_file,
                    focal_species=args.focal_species,
                    n_neighbors=args.n_neighbors,
                    session=session,
                    include_focal=True,
                    max_search=args.max_search,
                    min_sequences=args.min_sequences,
                )
                species_list = species_with_data
            else:
                neighbors = load_neighbors_from_tree(tree_file, args.focal_species, args.n_neighbors)
                species_list = [args.focal_species] + neighbors
        else:
            species_list = [args.focal_species]

        # Track actual number of focal clade species (with data)
        n_focal_clade = len(species_list)

        # Add outgroup samples - either by order (preferred) or by distance
        if args.outgroups_per_order > 0:
            # Order-based sampling (more thorough)
            if not args.require_ncbi:
                session = requests.Session()

            # Parse orders to sample
            if args.orders:
                orders_to_sample = [o.strip().upper() for o in args.orders.split(",")]
            else:
                # Default orders (excluding focal order, which is handled in the function)
                orders_to_sample = ["PRIMATES", "RODENTIA", "CARNIVORA", "ARTIODACTYLA", "CETACEA", "PERISSODACTYLA"]

            outgroups, order_counts = find_outgroups_by_order(
                tree_file=tree_file,
                focal_species=args.focal_species,
                n_per_order=args.outgroups_per_order,
                orders_to_sample=orders_to_sample,
                session=session,
                exclude_species=species_list,
                max_search_per_order=args.max_search,
                min_sequences=args.min_sequences,
            )
            species_list = species_list + outgroups

        elif args.outgroup_samples > 0:
            # Simple distance-based sampling
            if not args.require_ncbi:
                session = requests.Session()
            print(f"\n--- Adding {args.outgroup_samples} outgroup samples ---")
            outgroups, _ = find_outgroups_with_ncbi_data(
                tree_file=tree_file,
                focal_species=args.focal_species,
                n_outgroups=args.outgroup_samples,
                session=session,
                exclude_species=species_list,
                max_search=args.max_search,
                min_sequences=args.min_sequences,
            )
            print(f"Outgroups to test: {outgroups}")
            species_list = species_list + outgroups

    # Handle --target-species case (no focal clade tracking needed)
    if args.target_species:
        n_focal_clade = len(species_list)  # All are "focal" in single-target mode

    print(f"\n=== Final species list ({len(species_list)} species) ===")
    print(f"    Focal clade: {n_focal_clade} species")
    print(f"    Outgroups: {len(species_list) - n_focal_clade} species")
    for i, sp in enumerate(species_list):
        label = "(focal)" if i == 0 else "(neighbor)" if i < n_focal_clade else "(outgroup)"
        print(f"  {i+1}. {sp} {label}")

    results = validate_gene_multi(
        gene_symbol=gene_symbol,
        species_list=species_list,
        reference_source=args.reference_source,
        use_protein=args.protein,
        reciprocal_max_len=args.reciprocal_max_length,
    )

    # Add sample_type to results (focal, neighbor, or outgroup)
    for i, result in enumerate(results):
        if i == 0:
            result["sample_type"] = "focal"
        elif i < n_focal_clade:
            result["sample_type"] = "neighbor"
        else:
            result["sample_type"] = "outgroup"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.output:
        output_file = Path(args.output)
        # Create parent directory if it doesn't exist
        output_file.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Use focal species or target species for filename
        species_for_name = args.focal_species if args.focal_species else args.target_species
        species_clean = species_for_name.replace(" ", "_").lower()
        output_file = OUTPUT_DIR / f"{gene_symbol}_{species_clean}_validation.csv"

    fieldnames = [
        "gene_symbol",
        "target_species",
        "sample_type",
        "status",
        # Genome quality metrics
        "nuccore_count",
        "assembly_accession",
        "assembly_level",
        "scaffold_n50",
        "contig_n50",
        "total_length",
        # BLAST results
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

    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = {key: result.get(key, "") for key in fieldnames}
            if "reciprocal_details" in result:
                row["reciprocal_top_hit_accession"] = result["reciprocal_details"].get("top_hit_accession", "")
                row["reciprocal_top_hit_identity"] = result["reciprocal_details"].get("top_hit_identity", "")
                row["reciprocal_top_hit_evalue"] = result["reciprocal_details"].get("top_hit_evalue", "")
                row["reciprocal_top_hit_definition"] = result["reciprocal_details"].get("top_hit_definition", "")
            writer.writerow(row)

    print(f"Results saved to: {output_file}")

    if args.write_json:
        json_file = output_file.with_suffix(".json")
        with open(json_file, "w") as f:
            json.dump(results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
