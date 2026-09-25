# Reciprocal BLAST Gene Validation (API-only)

This tool validates gene gain/loss predictions using NCBI BLAST API calls and reciprocal best-hit logic.

## Overview

The validation workflow:
1. Fetches a human reference sequence from UniProt (protein) or Ensembl (CDS).
2. Runs BLAST (tblastn or blastn) against the target taxon via NCBI API.
3. Runs reciprocal BLAST of the top hit back to human.
4. Accepts homologs only if the original gene is the top reciprocal hit.

## Installation

```bash
pip install -r requirements.txt
```

Requirements:
- requests >= 2.31

## Usage

### Reciprocal BLAST Run

```bash
python validate_gene.py \
  --gene NKG7 \
  --target-species "Myotis davidii" \
  --reference-source uniprot \
  --protein
```

Run focal species plus N neighbors (from file):

```bash
python validate_gene.py \
  --gene NKG7 \
  --focal-species "Myotis davidii" \
  --n-neighbors 8 \
  --neighbors-file neighbors.txt \
  --reference-source uniprot \
  --protein
```

Run focal species plus N neighbors (from tree by branch length):

```bash
python validate_gene.py \
  --gene NKG7 \
  --focal-species "Myotis davidii" \
  --n-neighbors 8 \
  --tree-file path/to/tree.newick \
  --reference-source uniprot \
  --protein
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--gene` | Human gene symbol | Required |
| `--target-species` | Target species name | Optional |
| `--focal-species` | Focal species name | Optional |
| `--n-neighbors` | Number of neighbors to include | 0 |
| `--neighbors-file` | File with one neighbor per line | None |
| `--tree-file` | Newick tree file (requires biopython) | None |
| `--reference-source` | auto/uniprot/ensembl | auto |
| `--protein` | Use protein query (tblastn) | False |
| `--reciprocal-max-length` | Max reciprocal sequence length | 8000 |
| `--output` | Output CSV file | Auto-generated |
| `--write-json` | Write JSON debug output | False |

### Batch Validation

Create a TSV file with columns: `gene_symbol`, `focal_species`, `n_neighbors`

```bash
python batch_validate.py --input genes_to_validate.tsv --output-dir results/ --delay 60
```

Batch runs produce per-gene CSVs plus `batch_validation_summary.csv`.

Example input file (`example_genes.tsv`):
```
gene_symbol	focal_species	n_neighbors
RAB40B	Sciurus carolinensis	10
NKG7	Myotis davidii	8
```

## Output

Results are saved to the specified `--output` path and include identity scores
for forward and reciprocal hits plus reciprocal pass/fail.

### Example Output (CSV)

```
gene_symbol,target_species,status,reference_source,reference_id,reference_length,n_hits,ecnc,duplicated_call,top_hit_accession,top_hit_identity,top_hit_evalue,top_hit_definition,reciprocal_status,reciprocal_top_hit_accession,reciprocal_top_hit_identity,reciprocal_top_hit_evalue,reciprocal_top_hit_definition
NKG7,ENSG00000105374,Homo sapiens,Myotis davidii,Myotis_davidii,found,5,3,2,True,XP_12345,82.4,1e-50,Example hit,confirmed,NP_67890,91.2,1e-80,Human NK cell granule protein
```

## Interpreting Results

### Status Values

| Status | Meaning |
|--------|---------|
| `found` | Forward hit found (reciprocal status in `reciprocal_status`) |
| `no_hits` | Forward search found no hits |
| `no_genome` | Target genome missing or unavailable |

### ECNC duplication rule

ECNC is computed as $\text{ECNC} = \frac{\sum \text{aligned lengths}}{\text{query length}}$.

A gene is called **duplicated** when:
- total hit count $\ge 2$, **and**
- ECNC $\ge 1.5$ (Bowhead criteria)

### Duplication Detection

A gene is flagged as "likely_duplicated" when:
- Multiple distinct BLAST hits with E < 1e-10 and >70% identity
- Hits map to different genomic loci (different LOC IDs or genes)

**Note:** Multiple transcript variants of the same gene are NOT counted as duplications.

## Caveats

1. **Paralog confusion**: Reciprocal best hit helps, but close paralogs can still occur.

2. **No genome data**: Missing or incomplete target genomes yield `no_genome` or `no_hits` and should not be counted as true losses.

4. **Assembly quality**: Poor genome assemblies may miss genes that are actually present.

## Examples

### Test RAB40B Loss in Squirrels

TOGA predicts RAB40B is present only in *Sciurus carolinensis* and lost in other squirrels:

```bash
python validate_gene.py --gene RAB40B --focal-species "Sciurus carolinensis" -n 10
```

### Test NKG7 Duplication in Myotis Bats

TOGA predicts NKG7 is duplicated in some Myotis species:

```bash
python validate_gene.py --gene NKG7 --focal-species "Myotis davidii" -n 8 --check-duplication

# Use UniProt human reference with reciprocal BLAST back to both human and focal
python validate_gene.py --gene NKG7 --focal-species "Myotis davidii" -n 8 --reference-source uniprot --reciprocal --reciprocal-targets human,focal

# Run both reference methods (human + focal) and write separate CSVs
python validate_gene.py --gene NKG7 --focal-species "Myotis davidii" -n 8 --reference-method both
```

## Trait Correlation Prep (BayesTraits)

Use the helper script to merge TOGA gene gain/loss calls with trait data (e.g.,
`max_longevity_y`) and build BayesTraits input files. It writes:

- `merged_trait_data.csv`
- `bayestraits_input.txt`
- `bayestraits_commands.example.txt`

```bash
python toga_trait_analysis.py \
  --gene-states example_gene_states.tsv \
  --traits example_traits.tsv \
  --trait-column max_longevity_y \
  --discretize-trait
```

To run BayesTraits, provide the executable path, tree file, and a command file:

```bash
python toga_trait_analysis.py \
  --gene-states gene_states.tsv \
  --traits traits.tsv \
  --trait-column max_longevity_y \
  --tree /path/to/tree.newick \
  --bayestraits-path /path/to/BayesTraits \
  --bayestraits-commands bayestraits_commands.txt
```

## File Structure

```
scripts/copy_num_BLAT_check/
├── README.md                 # This file
├── validate_gene.py          # Main validation script
├── batch_validate.py         # Batch processing script
├── toga_trait_analysis.py     # BayesTraits prep script
├── run_validation.sh         # Shell wrapper
├── requirements.txt          # Python dependencies
└── example_genes.tsv         # Example batch input
├── example_gene_states.tsv    # Example gene-state input
├── example_traits.tsv         # Example trait input

output/copy_num_BLAT_check/
└── [gene]_[species]_validation.csv   # Results files
```
