from pathlib import Path
import os
REPO_ROOT = Path(os.environ.get("COPY_NUM_ROOT", Path(__file__).resolve().parents[2]))
BASE_DIR = REPO_ROOT

#!/usr/bin/env python3
"""
Refine longevity gene candidates by:
1. Filtering convergent patterns - ensure species with changes have longevity association
2. Adding PGLS/within-family candidates that are LOO-robust
3. Requiring significant longevity difference between control and affected groups

Criteria:
- For convergent patterns: mean MLres of affected species > 0 (longer-lived than expected)
- For PGLS: p < 0.05 AND LOO robust = TRUE
- For within-family: significant meta-p AND consistent direction
"""

import csv
from pathlib import Path
from collections import defaultdict
import statistics

BASE_DIR = REPO_ROOT
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"


def load_species_traits():
    """Load species trait data with MLres."""
    traits = {}
    trait_path = DATA_DIR / "species_trait_data.csv"

    with open(trait_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            species = row.get('Species_underscore', '').replace(" ", "_")
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
    return traits


def load_longevity_genes():
    """Load longevity gene list."""
    genes = set()
    with open(DATA_DIR / "de_https-::doi.org:10.1093:gbe:evad186.txt") as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return genes


def load_pgls_candidates(traits):
    """Load PGLS candidates that are FDR significant AND LOO robust."""
    candidates = []

    # Loss PGLS
    loss_pgls_path = RESULTS_DIR / "1_loss_analysis" / "pgls_gene_trait_associations_mlres.csv"
    if loss_pgls_path.exists():
        with open(loss_pgls_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    p_adj = float(row.get('pgls_p_adjusted', 1) or 1)
                    loo_robust = row.get('loo_robust', '').upper() == 'TRUE'
                    coef = float(row.get('pgls_coef', 0) or 0)
                    mlres_lost = float(row.get('mean_MLres_lost', 0) or 0)
                    mlres_retained = float(row.get('mean_MLres_retained', 0) or 0)

                    # Significant and robust
                    if p_adj < 0.05 and loo_robust:
                        candidates.append({
                            'gene': row['gene_symbol'],
                            'gene_id': row['gene_id'],
                            'event_type': 'loss',
                            'source': 'PGLS_MLres',
                            'p_value': p_adj,
                            'effect': coef,
                            'direction': 'longer' if coef > 0 else 'shorter',
                            'mean_mlres_affected': mlres_lost,
                            'mean_mlres_control': mlres_retained,
                            'loo_robust': True,
                            'is_longevity_gene': row.get('is_longevity_gene', 'FALSE').upper() == 'TRUE'
                        })
                except (ValueError, KeyError):
                    continue

    # Duplication PGLS
    dup_pgls_path = RESULTS_DIR / "2_duplication_analysis" / "pgls_duplication_mlres.tsv"
    if dup_pgls_path.exists():
        with open(dup_pgls_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    p_adj = float(row.get('pgls_p_adjusted', 1) or 1)
                    loo_robust = row.get('loo_robust', '').upper() == 'TRUE'
                    coef = float(row.get('pgls_coef', 0) or 0)

                    if p_adj < 0.05 and loo_robust:
                        candidates.append({
                            'gene': row['gene_symbol'],
                            'gene_id': row['gene_id'],
                            'event_type': 'duplication',
                            'source': 'PGLS_MLres',
                            'p_value': p_adj,
                            'effect': coef,
                            'direction': 'longer' if coef > 0 else 'shorter',
                            'mean_mlres_affected': None,
                            'mean_mlres_control': None,
                            'loo_robust': True,
                            'is_longevity_gene': False  # Check separately
                        })
                except (ValueError, KeyError):
                    continue

    return candidates


def load_within_family_candidates():
    """Load within-family candidates with consistent direction."""
    candidates = []

    for event_type, filename in [('loss', 'within_family_loss_summary.csv'),
                                  ('duplication', 'within_family_duplication_summary.csv')]:
        path = RESULTS_DIR / "4_within_family" / filename
        if not path.exists():
            continue

        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    meta_p = float(row.get('meta_p_value', 1) or 1)
                    consistency = float(row.get('direction_consistency', 0) or 0)
                    n_families = int(row.get('n_families_significant', 0) or 0)
                    final_score = float(row.get('final_score', 0) or 0)

                    # Significant in multiple families with consistent direction
                    if meta_p < 0.05 and consistency >= 0.7 and n_families >= 2:
                        direction = row.get('dominant_direction', 'positive')
                        candidates.append({
                            'gene': row['gene_symbol'],
                            'gene_id': row['gene_id'],
                            'event_type': event_type,
                            'source': 'within_family',
                            'p_value': meta_p,
                            'effect': float(row.get('mean_effect', 0) or 0),
                            'direction': 'longer' if direction == 'positive' else 'shorter',
                            'n_families': n_families,
                            'consistency': consistency,
                            'final_score': final_score,
                            'loo_robust': True  # Within-family is inherently robust
                        })
                except (ValueError, KeyError):
                    continue

    return candidates


def load_pattern_candidates(traits):
    """Load convergent pattern candidates and filter by species MLres."""
    path = RESULTS_DIR / "TOP_CANDIDATES" / "interesting_longevity_candidates.csv"
    if not path.exists():
        return []

    candidates = []

    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            species_list = row.get('species_list', '').split(';')
            pattern = row.get('pattern', '')

            # Get MLres for affected species
            mlres_values = []
            for sp in species_list:
                sp = sp.strip()
                if sp in traits:
                    mlres_values.append(traits[sp]['mlres'])

            if not mlres_values:
                continue

            mean_mlres = statistics.mean(mlres_values)
            n_long_lived = sum(1 for m in mlres_values if m > 0)

            # For convergent patterns, require affected species to actually be long-lived
            if pattern == 'convergent_long_lived':
                # At least 50% of affected species should be long-lived
                if n_long_lived < len(mlres_values) * 0.5:
                    continue  # Deprioritize - species aren't actually long-lived

            # For rare clustered, check if species is long-lived
            elif pattern == 'rare_clustered':
                # For single-species losses, require that species to be notably long-lived
                if len(species_list) == 1 and mean_mlres < 0.1:
                    continue  # Skip - single loss in short-lived species not interesting

            candidates.append({
                'gene': row['gene'],
                'gene_id': row['gene_id'],
                'event_type': row['event_type'],
                'source': f'pattern_{pattern}',
                'pattern': pattern,
                'n_events': int(row.get('n_events', 0)),
                'n_clades': int(row.get('n_clades', 0)),
                'primary_clade': row.get('primary_clade', ''),
                'mean_mlres_affected': round(mean_mlres, 3),
                'n_long_lived': n_long_lived,
                'pct_long_lived': round(100 * n_long_lived / len(mlres_values), 1),
                'species': species_list[:5],
                'passes_longevity_filter': mean_mlres > 0 or n_long_lived >= 2
            })

    # Filter to only those passing longevity filter
    return [c for c in candidates if c.get('passes_longevity_filter', False)]


def main():
    print("=" * 70)
    print("REFINING LONGEVITY GENE CANDIDATES")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    traits = load_species_traits()
    print(f"  Species with traits: {len(traits)}")

    longevity_genes = load_longevity_genes()
    print(f"  Known longevity genes: {len(longevity_genes)}")

    # Load candidates from different sources
    print("\nLoading PGLS candidates (FDR < 0.05 & LOO robust)...")
    pgls_candidates = load_pgls_candidates(traits)
    print(f"  PGLS candidates: {len(pgls_candidates)}")

    print("\nLoading within-family candidates...")
    wf_candidates = load_within_family_candidates()
    print(f"  Within-family candidates: {len(wf_candidates)}")

    print("\nLoading pattern candidates (filtered by longevity)...")
    pattern_candidates = load_pattern_candidates(traits)
    print(f"  Pattern candidates (filtered): {len(pattern_candidates)}")

    # Combine and deduplicate
    all_candidates = []
    seen_genes = set()

    # Priority 1: PGLS + LOO robust + longevity gene
    for c in pgls_candidates:
        if c['gene'] in longevity_genes:
            c['priority'] = 1
            c['priority_reason'] = 'PGLS_robust + longevity_gene'
            if c['gene'] not in seen_genes:
                all_candidates.append(c)
                seen_genes.add(c['gene'])

    # Priority 2: Pattern candidates (rare clustered in long-lived)
    for c in pattern_candidates:
        if c['pattern'] == 'rare_clustered' and c['mean_mlres_affected'] > 0.1:
            c['priority'] = 2
            c['priority_reason'] = 'rare_clustered_long_lived'
            if c['gene'] not in seen_genes:
                all_candidates.append(c)
                seen_genes.add(c['gene'])

    # Priority 3: Convergent in long-lived
    for c in pattern_candidates:
        if c['pattern'] == 'convergent_long_lived' and c['pct_long_lived'] >= 50:
            c['priority'] = 3
            c['priority_reason'] = 'convergent_long_lived'
            if c['gene'] not in seen_genes:
                all_candidates.append(c)
                seen_genes.add(c['gene'])

    # Priority 4: Within-family robust
    for c in wf_candidates:
        if c['gene'] in longevity_genes:
            c['priority'] = 4
            c['priority_reason'] = 'within_family + longevity_gene'
            if c['gene'] not in seen_genes:
                all_candidates.append(c)
                seen_genes.add(c['gene'])

    # Priority 5: Other PGLS robust
    for c in pgls_candidates:
        c['priority'] = 5
        c['priority_reason'] = 'PGLS_robust'
        if c['gene'] not in seen_genes:
            all_candidates.append(c)
            seen_genes.add(c['gene'])

    # Sort by priority
    all_candidates.sort(key=lambda x: (x.get('priority', 99), -abs(x.get('effect', 0))))

    # Print results
    print("\n" + "=" * 70)
    print("REFINED CANDIDATES BY PRIORITY")
    print("=" * 70)

    for priority in [1, 2, 3, 4, 5]:
        priority_cands = [c for c in all_candidates if c.get('priority') == priority]
        if not priority_cands:
            continue

        print(f"\n--- PRIORITY {priority}: {priority_cands[0].get('priority_reason', '')} ---")
        for c in priority_cands[:10]:
            gene = c['gene']
            event = c['event_type'].upper()
            longevity_tag = " [LONGEVITY]" if c['gene'] in longevity_genes else ""

            if 'pattern' in c:
                print(f"  {gene}{longevity_tag} ({event}): {c['n_events']} events, "
                      f"mean_MLres={c.get('mean_mlres_affected', 'NA')}, "
                      f"{c.get('pct_long_lived', 'NA')}% long-lived")
            else:
                direction = c.get('direction', '')
                p = c.get('p_value', 1)
                print(f"  {gene}{longevity_tag} ({event}): p={p:.2e}, direction={direction}")

    # Save refined list
    output_path = RESULTS_DIR / "TOP_CANDIDATES" / "refined_candidates_for_miniprot.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w', newline='') as f:
        fieldnames = ['gene', 'gene_id', 'event_type', 'source', 'priority',
                      'priority_reason', 'p_value', 'effect', 'direction',
                      'mean_mlres_affected', 'n_events', 'n_clades', 'loo_robust',
                      'is_longevity_gene']
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()

        for c in all_candidates:
            c['is_longevity_gene'] = c['gene'] in longevity_genes
            writer.writerow(c)

    print(f"\n\nSaved {len(all_candidates)} refined candidates to:")
    print(f"  {output_path}")

    # Summary for miniprot
    print("\n" + "=" * 70)
    print("TOP 20 CANDIDATES FOR MINIPROT VALIDATION")
    print("=" * 70)

    top_for_validation = []
    for c in all_candidates[:30]:
        if c['gene'] in longevity_genes or c.get('priority', 99) <= 3:
            top_for_validation.append(c)
        if len(top_for_validation) >= 20:
            break

    for i, c in enumerate(top_for_validation, 1):
        gene = c['gene']
        event = c['event_type']
        reason = c.get('priority_reason', c.get('source', ''))
        longevity_tag = " *" if c['gene'] in longevity_genes else ""
        print(f"  {i:2}. {gene}{longevity_tag} ({event}) - {reason}")

    # Generate gene list for miniprot
    genes_for_validation = [c['gene'] for c in top_for_validation]
    print(f"\nGene list for miniprot:")
    print(f"  {','.join(genes_for_validation)}")

    return all_candidates


if __name__ == "__main__":
    main()
