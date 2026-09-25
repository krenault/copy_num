#!/bin/bash
# Minimal miniprot validation using ONLY already-cached genomes
# No new genome downloads - stays within disk space limit

set -e
cd "$(dirname "$0")"

GENOME_DIR="genomes"
RESULTS_DIR="results"
MIN_IDENTITY=80

mkdir -p "$RESULTS_DIR"

echo "============================================================"
echo "MINIMAL MINIPROT VALIDATION (CACHED GENOMES ONLY)"
echo "============================================================"
echo ""
echo "Using only pre-cached genomes to stay within disk space limit."
echo ""

# Species with already-downloaded genomes
CACHED_SPECIES=(
    "Homo_sapiens"           # Reference
    "Pan_troglodytes"        # Primate control
    "Macaca_mulatta"         # Primate control
    "Mus_musculus"           # Rodent control
    "Rattus_norvegicus"      # Rodent control
    "Sciurus_carolinensis"   # AFFECTED: ZMPSTE24 + HTRA2 loss
    "Myotis_myotis"          # Bat control
    "Pteropus_vampyrus"      # Bat control
    "Panthera_tigris"        # Carnivore control
    "Canis_lupus_familiaris" # Carnivore control (may be cached as Canis_lupus_baileyi)
)

SPECIES_STR=$(IFS=,; echo "${CACHED_SPECIES[*]}")

echo "Cached species (${#CACHED_SPECIES[@]}):"
printf "  %s\n" "${CACHED_SPECIES[@]}"
echo ""

# TOP CANDIDATES where we have affected OR control species cached

echo "============================================================"
echo "PRIORITY 1: RARE LOSSES WITH CACHED AFFECTED SPECIES"
echo "============================================================"
echo ""
echo "ZMPSTE24 - Loss in Sciurus_carolinensis (squirrel)"
echo "  Known PROGERIA gene - loss in long-lived squirrel very interesting!"
echo "HTRA2 - Loss in Sciurus_carolinensis"
echo "  Mitochondrial serine protease, Parkinson's association"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "ZMPSTE24,HTRA2" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_squirrel_losses.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes \
    --skip-download 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 2: CONVERGENT LOSSES IN LONG-LIVED (control check)"
echo "============================================================"
echo ""
echo "APOE - Losses across primates/bats (Alzheimer's gene)"
echo "GSTP1 - Detoxification/longevity gene"
echo "ATM - DNA damage response"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "APOE,GSTP1,ATM" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_losses.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes \
    --skip-download 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 3: CONVERGENT DUPLICATIONS (control check)"
echo "============================================================"
echo ""
echo "MTOR - Master metabolic regulator"
echo "CDKN1A - Cell cycle (p21)"
echo "APEX1 - DNA repair"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "MTOR,CDKN1A,APEX1" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_dups.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes \
    --skip-download 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 4: RARE CLUSTERED LOSSES (control species)"
echo "============================================================"
echo ""
echo "CDC42 - Loss in Myotis_brandtii (NOT cached - checking controls)"
echo "PDGFRB - Cell signaling"
echo "TXN - Thioredoxin antioxidant"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "CDC42,PDGFRB,TXN" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_rare_losses_control.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes \
    --skip-download 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "VALIDATION COMPLETE"
echo "============================================================"
echo ""
echo "Results:"
ls -la "$RESULTS_DIR"/miniprot_*.csv 2>/dev/null || echo "  Checking..."
echo ""
echo "Key findings to look for:"
echo "  1. ZMPSTE24 in Sciurus_carolinensis: should show 0 loci (TRUE loss)"
echo "  2. HTRA2 in Sciurus_carolinensis: should show 0 loci (TRUE loss)"
echo "  3. Control species should show 1 locus for most genes"
echo ""
echo "NOTE: To validate additional affected species (bats with CDC42 loss, etc.),"
echo "you would need to clear space and download those specific genomes."
