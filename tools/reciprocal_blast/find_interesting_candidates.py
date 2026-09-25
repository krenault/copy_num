from pathlib import Path
import os
REPO_ROOT = Path(os.environ.get("COPY_NUM_ROOT", Path(__file__).resolve().parents[2]))
BASE_DIR = REPO_ROOT

#!/usr/bin/env python3
"""
Find longevity genes with interesting loss/duplication patterns for miniprot validation.

Criteria for "interesting":
1. RARE CLUSTERED CHANGES (like CGAS in armadillos):
   - 1-5 losses/duplications in closely related species (same order)
   - High clade concentration

2. CONVERGENT CHANGES IN DIVERGED LONG-LIVED SPECIES:
   - Changes in species from different orders
   - Both have high longevity (MLres > 0)

3. EXCLUDE:
   - >15 changes all in one clade (too many, likely artifacts)
   - Genes with changes spread randomly across phylogeny
"""

import csv
from pathlib import Path
from collections import defaultdict

# Paths
BASE_DIR = REPO_ROOT
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"

def load_longevity_genes():
    """Load longevity gene list."""
    genes = set()
    with open(DATA_DIR / "de_https-::doi.org:10.1093:gbe:evad186.txt") as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return genes

def load_loss_matrix():
    """Load gene loss from presence/absence matrix (0 = lost, 1 = present)."""
    matrix = {}
    # Use presence/absence matrix - 0 means gene is LOST
    matrix_path = DATA_DIR / "loss_analysis" / "All_Species_Gene_PresenceAbsence.tsv"
    if not matrix_path.exists():
        print(f"  WARNING: Loss matrix not found at {matrix_path}")
        return {}, []

    with open(matrix_path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        species_cols = header[2:]  # Skip gene_id, gene_symbol

        for row in reader:
            if len(row) < 3:
                continue
            gene_id = row[0]
            gene_symbol = row[1]
            losses = {}
            for i, sp in enumerate(species_cols):
                sp_clean = sp.replace(".", "_")
                try:
                    # In presence/absence: 0 = LOST, 1 = present
                    # Convert to loss matrix: 1 = LOST, 0 = present
                    presence = int(row[i + 2])
                    losses[sp_clean] = 1 if presence == 0 else 0
                except (ValueError, IndexError):
                    losses[sp_clean] = 0
            matrix[gene_symbol] = {'id': gene_id, 'losses': losses}

    return matrix, species_cols

def load_duplication_matrix():
    """Load gene duplication binary matrix."""
    matrix = {}
    matrix_path = DATA_DIR / "duplication_dollo" / "All_Species_Gene_Duplication_Binary.tsv"
    if not matrix_path.exists():
        return {}, []

    with open(matrix_path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        species_cols = header[2:]

        for row in reader:
            if len(row) < 3:
                continue
            gene_id = row[0]
            gene_symbol = row[1]
            values = {}
            for i, sp in enumerate(species_cols):
                sp_clean = sp.replace(".", "_")
                try:
                    values[sp_clean] = int(row[i + 2])
                except (ValueError, IndexError):
                    values[sp_clean] = 0
            matrix[gene_symbol] = {'id': gene_id, 'duplications': values}

    return matrix, species_cols

def load_species_traits():
    """Load species trait data including order and MLres."""
    traits = {}
    trait_path = DATA_DIR / "species_trait_data.csv"

    if trait_path.exists():
        with open(trait_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Get species name - use Species_underscore column
                species = row.get('Species_underscore', row.get('Species', '')).replace(" ", "_")
                if species:
                    try:
                        mlres = float(row.get('MLres', 0) or 0)
                    except:
                        mlres = 0
                    try:
                        max_long = float(row.get('maximum_longevity_y', 0) or 0)
                    except:
                        max_long = 0

                    traits[species] = {
                        'order': row.get('order', 'Unknown'),
                        'mlres': mlres,
                        'max_longevity': max_long
                    }
        print(f"  Loaded traits for {len(traits)} species")
    else:
        print(f"  WARNING: Trait file not found at {trait_path}")

    return traits

def extract_order(species_name):
    """Extract order from species name like 'Myotis_myotis_VESPERTILIONIDAE_CHIROPTERA'."""
    parts = species_name.split("_")
    if len(parts) >= 4:
        # Order is typically last part
        potential_order = parts[-1].upper()
        if potential_order in ['PRIMATES', 'RODENTIA', 'CHIROPTERA', 'CARNIVORA',
                               'ARTIODACTYLA', 'CETACEA', 'PERISSODACTYLA', 'LAGOMORPHA',
                               'EULIPOTYPHLA', 'PHOLIDOTA', 'CINGULATA', 'PILOSA']:
            return potential_order
    return "Unknown"

def analyze_loss_patterns(loss_matrix, longevity_genes, species_traits):
    """Find longevity genes with interesting loss patterns."""
    candidates = []

    for gene_symbol, data in loss_matrix.items():
        if gene_symbol not in longevity_genes:
            continue

        losses = data['losses']
        species_with_loss = [sp for sp, val in losses.items() if val == 1]
        n_losses = len(species_with_loss)

        if n_losses == 0 or n_losses > 15:
            continue  # Skip genes with no losses or too many

        # Analyze clade distribution
        orders = defaultdict(list)
        for sp in species_with_loss:
            order = species_traits.get(sp, {}).get('order', extract_order(sp))
            orders[order].append(sp)

        # Calculate clade concentration
        max_in_one_clade = max(len(sps) for sps in orders.values()) if orders else 0
        n_clades = len(orders)

        # Pattern 1: Rare clustered (like CGAS) - most losses in one clade, few total
        if n_losses <= 5 and max_in_one_clade >= n_losses * 0.8:
            pattern = "rare_clustered"
            primary_clade = max(orders.keys(), key=lambda x: len(orders[x]))
            candidates.append({
                'gene': gene_symbol,
                'gene_id': data['id'],
                'event_type': 'loss',
                'pattern': pattern,
                'n_events': n_losses,
                'n_clades': n_clades,
                'primary_clade': primary_clade,
                'species': species_with_loss,
                'clade_concentration': max_in_one_clade / n_losses if n_losses > 0 else 0
            })

        # Pattern 2: Convergent in diverged long-lived (losses in different orders, high MLres)
        elif n_clades >= 2 and n_losses <= 10:
            long_lived_losses = []
            for sp in species_with_loss:
                mlres = species_traits.get(sp, {}).get('mlres', 0)
                if mlres > 0.1:  # Positive MLres = longer-lived than expected
                    long_lived_losses.append(sp)

            if len(long_lived_losses) >= 2:
                # Check if they're from different clades
                ll_orders = set(species_traits.get(sp, {}).get('order', extract_order(sp))
                               for sp in long_lived_losses)
                if len(ll_orders) >= 2:
                    pattern = "convergent_long_lived"
                    candidates.append({
                        'gene': gene_symbol,
                        'gene_id': data['id'],
                        'event_type': 'loss',
                        'pattern': pattern,
                        'n_events': n_losses,
                        'n_clades': n_clades,
                        'primary_clade': 'multiple',
                        'species': species_with_loss,
                        'long_lived_species': long_lived_losses,
                        'clade_concentration': max_in_one_clade / n_losses if n_losses > 0 else 0
                    })

    return candidates

def analyze_duplication_patterns(dup_matrix, longevity_genes, species_traits):
    """Find longevity genes with interesting duplication patterns."""
    candidates = []

    for gene_symbol, data in dup_matrix.items():
        if gene_symbol not in longevity_genes:
            continue

        dups = data['duplications']
        species_with_dup = [sp for sp, val in dups.items() if val == 1]
        n_dups = len(species_with_dup)

        if n_dups == 0:
            continue

        # Analyze clade distribution
        orders = defaultdict(list)
        for sp in species_with_dup:
            order = species_traits.get(sp, {}).get('order', extract_order(sp))
            orders[order].append(sp)

        max_in_one_clade = max(len(sps) for sps in orders.values()) if orders else 0
        n_clades = len(orders)

        # Skip if >15 in one clade (too many, not informative)
        if max_in_one_clade > 15:
            continue

        # Pattern 1: Rare clustered (1-5 dups mostly in one clade)
        if n_dups <= 5 and max_in_one_clade >= n_dups * 0.8:
            pattern = "rare_clustered"
            primary_clade = max(orders.keys(), key=lambda x: len(orders[x]))
            candidates.append({
                'gene': gene_symbol,
                'gene_id': data['id'],
                'event_type': 'duplication',
                'pattern': pattern,
                'n_events': n_dups,
                'n_clades': n_clades,
                'primary_clade': primary_clade,
                'species': species_with_dup,
                'clade_concentration': max_in_one_clade / n_dups if n_dups > 0 else 0
            })

        # Pattern 2: Convergent in diverged long-lived
        elif n_clades >= 2 and n_dups <= 15:
            long_lived_dups = []
            for sp in species_with_dup:
                mlres = species_traits.get(sp, {}).get('mlres', 0)
                if mlres > 0.1:
                    long_lived_dups.append(sp)

            if len(long_lived_dups) >= 2:
                ll_orders = set(species_traits.get(sp, {}).get('order', extract_order(sp))
                               for sp in long_lived_dups)
                if len(ll_orders) >= 2:
                    pattern = "convergent_long_lived"
                    candidates.append({
                        'gene': gene_symbol,
                        'gene_id': data['id'],
                        'event_type': 'duplication',
                        'pattern': pattern,
                        'n_events': n_dups,
                        'n_clades': n_clades,
                        'primary_clade': 'multiple',
                        'species': species_with_dup,
                        'long_lived_species': long_lived_dups,
                        'clade_concentration': max_in_one_clade / n_dups if n_dups > 0 else 0
                    })

    return candidates

def main():
    print("=" * 70)
    print("FINDING INTERESTING LONGEVITY GENE CANDIDATES FOR MINIPROT VALIDATION")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    longevity_genes = load_longevity_genes()
    print(f"  Longevity genes: {len(longevity_genes)}")

    loss_matrix, loss_species = load_loss_matrix()
    print(f"  Loss matrix: {len(loss_matrix)} genes")

    dup_matrix, dup_species = load_duplication_matrix()
    print(f"  Duplication matrix: {len(dup_matrix)} genes")

    species_traits = load_species_traits()
    print(f"  Species with traits: {len(species_traits)}")

    # Analyze patterns
    print("\nAnalyzing loss patterns...")
    loss_candidates = analyze_loss_patterns(loss_matrix, longevity_genes, species_traits)
    print(f"  Found {len(loss_candidates)} interesting loss candidates")

    print("\nAnalyzing duplication patterns...")
    dup_candidates = analyze_duplication_patterns(dup_matrix, longevity_genes, species_traits)
    print(f"  Found {len(dup_candidates)} interesting duplication candidates")

    # Combine and sort
    all_candidates = loss_candidates + dup_candidates

    # Sort by pattern (rare_clustered first), then by number of events
    all_candidates.sort(key=lambda x: (
        0 if x['pattern'] == 'rare_clustered' else 1,
        x['n_events']
    ))

    # Print results
    print("\n" + "=" * 70)
    print("RARE CLUSTERED CHANGES (like CGAS)")
    print("=" * 70)

    rare_clustered = [c for c in all_candidates if c['pattern'] == 'rare_clustered']
    for c in rare_clustered[:20]:
        print(f"\n{c['gene']} ({c['event_type'].upper()})")
        print(f"  Events: {c['n_events']} | Clade: {c['primary_clade']}")
        print(f"  Species: {', '.join(c['species'][:5])}")

    print("\n" + "=" * 70)
    print("CONVERGENT CHANGES IN DIVERGED LONG-LIVED SPECIES")
    print("=" * 70)

    convergent = [c for c in all_candidates if c['pattern'] == 'convergent_long_lived']
    for c in convergent[:20]:
        print(f"\n{c['gene']} ({c['event_type'].upper()})")
        print(f"  Events: {c['n_events']} | Clades: {c['n_clades']}")
        print(f"  Long-lived species: {', '.join(c.get('long_lived_species', [])[:5])}")

    # Save to file
    output_path = RESULTS_DIR / "TOP_CANDIDATES" / "interesting_longevity_candidates.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        fieldnames = ['gene', 'gene_id', 'event_type', 'pattern', 'n_events',
                     'n_clades', 'primary_clade', 'clade_concentration',
                     'reason_flagged', 'species_list', 'long_lived_species_list']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for c in all_candidates:
            c['species_list'] = ';'.join(c.get('species', [])[:10])
            c['long_lived_species_list'] = ';'.join(c.get('long_lived_species', [])[:10])

            # Generate human-readable reason
            if c['pattern'] == 'rare_clustered':
                c['reason_flagged'] = (
                    f"Rare clustered {c['event_type']}: {c['n_events']} events "
                    f"concentrated in {c['primary_clade']} "
                    f"({c['clade_concentration']:.0%} in one clade)"
                )
            elif c['pattern'] == 'convergent_long_lived':
                ll_species = c.get('long_lived_species', [])[:3]
                c['reason_flagged'] = (
                    f"Convergent {c['event_type']} in {c['n_clades']} diverged clades; "
                    f"{len(c.get('long_lived_species', []))} long-lived species affected "
                    f"({', '.join(ll_species)}{'...' if len(c.get('long_lived_species', [])) > 3 else ''})"
                )
            else:
                c['reason_flagged'] = f"{c['pattern']}: {c['n_events']} {c['event_type']} events"

            writer.writerow(c)

    print(f"\n\nResults saved to: {output_path}")
    print(f"Total candidates: {len(all_candidates)}")

    # Print top candidates for miniprot
    print("\n" + "=" * 70)
    print("TOP CANDIDATES FOR MINIPROT VALIDATION")
    print("=" * 70)

    # Select best candidates
    top_for_miniprot = []

    # Add top rare_clustered
    for c in rare_clustered[:5]:
        top_for_miniprot.append(c)

    # Add top convergent
    for c in convergent[:5]:
        top_for_miniprot.append(c)

    # Check if CGAS is in there (the reference example)
    cgas_found = any(c['gene'] == 'CGAS' for c in all_candidates)
    if cgas_found:
        print("\n✓ CGAS found in candidates (reference example)")
    else:
        print("\n✗ CGAS not found - checking raw data...")

    print("\nGenes for miniprot validation:")
    for c in top_for_miniprot:
        print(f"  - {c['gene']}: {c['event_type']}, {c['n_events']} events in {c.get('primary_clade', 'multiple clades')}")

    return all_candidates


if __name__ == "__main__":
    main()
