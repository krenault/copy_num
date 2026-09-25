#!/bin/bash
# Overnight batch miniprot validation
# Validates each candidate gene separately, outputs to individual CSVs
#
# Usage: ./run_overnight_validation.sh [candidate_genes_file]
# Default: uses CANDIDATE_GENES array below

set -e
cd "$(dirname "$0")"

GENOME_DIR="genomes"
RESULTS_DIR="results/overnight_$(date +%Y%m%d)"
MIN_IDENTITY=80
LOG_FILE="$RESULTS_DIR/validation.log"

mkdir -p "$RESULTS_DIR"

# Default candidate genes (edit this list or pass a file)
# Template for good candidates: CLU in sloths - rare loss in specific clade
CANDIDATE_GENES=(
    # PRIORITY: CLU in sloths (excellent example - rare clade-specific loss)
    "CLU"
    # Loss candidates (squirrel)
    "ZMPSTE24"
    "HTRA2"
    # Convergent losses
    "APOE"
    "GSTP1"
    "ATM"
    # cGAS (armadillo)
    "CGAS"
    # Duplications to validate
    "MTOR"
    "CDKN1A"
    "APEX1"
    # Additional candidates
    "CDC42"
    "PDGFRB"
    "TXN"
)

# Species with cached genomes
CACHED_SPECIES=(
    "Homo_sapiens"
    "Pan_troglodytes"
    "Macaca_mulatta"
    "Mus_musculus"
    "Rattus_norvegicus"
    "Sciurus_carolinensis"
    "Myotis_myotis"
    "Pteropus_vampyrus"
    "Panthera_tigris"
    "Canis_lupus_familiaris"
)

# Additional species to download (comment out to skip)
DOWNLOAD_SPECIES=(
    "Dasypus_novemcinctus"  # Armadillo for CGAS
    "Choloepus_hoffmanni"   # Two-toed sloth for CLU
    "Bradypus_torquatus"    # Three-toed sloth for CLU
    # "Tolypeutes_matacus"    # Three-banded armadillo
)

# Read genes from file if provided
if [[ -n "$1" ]] && [[ -f "$1" ]]; then
    readarray -t CANDIDATE_GENES < "$1"
    echo "Loaded ${#CANDIDATE_GENES[@]} genes from $1"
fi

ALL_SPECIES=("${CACHED_SPECIES[@]}" "${DOWNLOAD_SPECIES[@]}")
SPECIES_STR=$(IFS=,; echo "${ALL_SPECIES[*]}")

echo "============================================================" | tee "$LOG_FILE"
echo "OVERNIGHT MINIPROT VALIDATION" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Validating ${#CANDIDATE_GENES[@]} genes:" | tee -a "$LOG_FILE"
printf "  %s\n" "${CANDIDATE_GENES[@]}" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Species (${#ALL_SPECIES[@]}):" | tee -a "$LOG_FILE"
printf "  %s\n" "${ALL_SPECIES[@]}" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Results directory: $RESULTS_DIR" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Process each gene individually
for GENE in "${CANDIDATE_GENES[@]}"; do
    echo "============================================================" | tee -a "$LOG_FILE"
    echo "Processing: $GENE" | tee -a "$LOG_FILE"
    echo "Time: $(date +%H:%M:%S)" | tee -a "$LOG_FILE"
    echo "============================================================" | tee -a "$LOG_FILE"

    OUTPUT_CSV="$RESULTS_DIR/${GENE}_validation.csv"

    python3 validate_candidates_miniprot.py \
        --genes "$GENE" \
        --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
        --species "$SPECIES_STR" \
        --output "$OUTPUT_CSV" \
        --genome-dir "$GENOME_DIR" \
        --min-n50 5000000 \
        --min-identity "$MIN_IDENTITY" \
        --keep-genomes 2>&1 | tee -a "$LOG_FILE" || echo "  ⚠ $GENE validation failed" | tee -a "$LOG_FILE"

    if [[ -f "$OUTPUT_CSV" ]]; then
        echo "  ✓ Saved: $OUTPUT_CSV" | tee -a "$LOG_FILE"
    fi
    echo "" | tee -a "$LOG_FILE"
done

echo "============================================================" | tee -a "$LOG_FILE"
echo "VALIDATION COMPLETE" | tee -a "$LOG_FILE"
echo "Finished: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Results:" | tee -a "$LOG_FILE"
ls -la "$RESULTS_DIR"/*.csv 2>/dev/null | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Storage after validation:" | tee -a "$LOG_FILE"
df -h . | tail -1 | tee -a "$LOG_FILE"
