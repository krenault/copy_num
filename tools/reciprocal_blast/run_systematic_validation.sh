#!/bin/bash
# Systematic miniprot validation of top TOGA duplication candidates
# Validates across multiple mammalian clades with good genome quality

set -e
cd "$(dirname "$0")"

# Configuration
GENOME_DIR="genomes"
RESULTS_DIR="results"
MIN_IDENTITY=80  # Strict ortholog threshold
MIN_N50=10000000  # 10Mb scaffold N50 minimum for quality

mkdir -p "$GENOME_DIR" "$RESULTS_DIR"

echo "============================================================"
echo "SYSTEMATIC MINIPROT VALIDATION OF TOGA DUPLICATION CANDIDATES"
echo "============================================================"
echo ""
echo "Identity threshold: ${MIN_IDENTITY}%"
echo "Minimum scaffold N50: ${MIN_N50}"
echo ""

# TOP CANDIDATE GENES FOR VALIDATION
# Selected based on:
# 1. High interest scores from duplication analysis
# 2. Clade-specific patterns (Chiroptera, Primates)
# 3. Longevity gene associations
# 4. Prior validation attempts (MINDY3, NKG7, TTC36)

GENES=(
    # Longevity-associated genes with duplication signals
    "POT1"      # Telomere protection, interest_score=75
    "XRCC1"     # DNA repair, interest_score=75
    "CISD2"     # Mitochondrial iron-sulfur, longevity gene

    # Chiroptera-concentrated duplications
    "HIGD2A"    # Hypoxia response, Chiroptera clade_concentration=13
    "TSGA10"    # Testis-specific, Chiroptera clade_concentration=7.29

    # Clade-specific from prior analysis
    "MINDY3"    # Deubiquitinase, Chiroptera-specific per TOGA
    "USP32"     # Ubiquitin protease, Primate-specific
    "PPP1R2B"   # Phosphatase regulator, Primate-specific

    # Primate-concentrated with high duplication
    "PPP2CA"    # Protein phosphatase, 78 species, Primate clade_concentration=78
    "VDAC1"     # Mitochondrial channel, 72 species, Primate
)

# SPECIES SELECTION - Broad clade coverage with high-quality genomes
# Format: "Species_name:Order"
SPECIES=(
    # Chiroptera (bats) - focal clade for many candidates
    "Myotis_myotis:Chiroptera"
    "Myotis_lucifugus:Chiroptera"
    "Pteropus_vampyrus:Chiroptera"
    "Rhinolophus_ferrumequinum:Chiroptera"

    # Primates - outgroup with many duplications
    "Homo_sapiens:Primates"
    "Pan_troglodytes:Primates"
    "Macaca_mulatta:Primates"

    # Rodentia - short-lived reference
    "Mus_musculus:Rodentia"
    "Rattus_norvegicus:Rodentia"

    # Carnivora - long-lived outgroup
    "Panthera_tigris:Carnivora"
    "Canis_lupus_familiaris:Carnivora"

    # Artiodactyla/Cetacea - diverse lifespans
    "Bos_taurus:Artiodactyla"
    "Sus_scrofa:Artiodactyla"

    # Perissodactyla - long-lived
    "Equus_caballus:Perissodactyla"
)

# Extract just species names for command
SPECIES_NAMES=""
for sp in "${SPECIES[@]}"; do
    name="${sp%%:*}"
    if [ -z "$SPECIES_NAMES" ]; then
        SPECIES_NAMES="$name"
    else
        SPECIES_NAMES="$SPECIES_NAMES,$name"
    fi
done

echo "Genes to validate: ${#GENES[@]}"
for g in "${GENES[@]}"; do
    echo "  - $g"
done
echo ""
echo "Species (${#SPECIES[@]} total):"
for sp in "${SPECIES[@]}"; do
    echo "  - ${sp%%:*} (${sp##*:})"
done
echo ""

# Run validation for each gene
for GENE in "${GENES[@]}"; do
    echo ""
    echo "============================================================"
    echo "Validating: $GENE"
    echo "============================================================"

    OUTPUT_FILE="${RESULTS_DIR}/miniprot_${GENE}.csv"

    python3 validate_candidates_miniprot.py \
        --genes "$GENE" \
        --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
        --species "$SPECIES_NAMES" \
        --output "$OUTPUT_FILE" \
        --genome-dir "$GENOME_DIR" \
        --min-n50 "$MIN_N50" \
        --min-identity "$MIN_IDENTITY" \
        --keep-genomes

    echo "Results saved to: $OUTPUT_FILE"
done

echo ""
echo "============================================================"
echo "VALIDATION COMPLETE"
echo "============================================================"
echo ""
echo "Results files:"
ls -la "$RESULTS_DIR"/miniprot_*.csv 2>/dev/null || echo "  (no results yet)"
echo ""
echo "To view summary:"
echo "  cat ${RESULTS_DIR}/miniprot_*.csv | column -t -s,"
