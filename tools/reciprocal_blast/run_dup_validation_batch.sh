#!/usr/bin/env zsh
# Batch miniprot validation for DUPLICATION candidates
# Outputs separate CSV for each gene
#
# Usage: ./run_dup_validation_batch.sh <gene_list_file> [max_genes]
# Example: ./run_dup_validation_batch.sh top25_dup_candidates.txt

set -e
cd "$(dirname "$0")"

GENE_FILE="${1:-top25_dup_candidates.txt}"
MAX_GENES="${2:-0}"  # 0 = all genes

GENOME_DIR="genomes"
RESULTS_DIR="results/dup_validation_$(date +%Y%m%d_%H%M)"
MIN_IDENTITY=70  # Slightly lower for detecting diverged duplicates
LOG_FILE="$RESULTS_DIR/validation.log"

mkdir -p "$RESULTS_DIR"

# Species to validate
ALL_SPECIES=(
    # Cached species
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
    # Sloths
    "Choloepus_didactylus"
    # Armadillo
    "Dasypus_novemcinctus"
)

SPECIES_STR=$(IFS=,; echo "${ALL_SPECIES[*]}")

# Read genes from file
if [[ ! -f "$GENE_FILE" ]]; then
    echo "Error: Gene file not found: $GENE_FILE"
    exit 1
fi

# Read genes into array (portable)
ALL_GENES=()
while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -n "$line" ]] && ALL_GENES+=("$line")
done < "$GENE_FILE"
TOTAL_GENES=${#ALL_GENES[@]}

if [[ $MAX_GENES -gt 0 ]] && [[ $MAX_GENES -lt $TOTAL_GENES ]]; then
    GENES=("${ALL_GENES[@]:0:$MAX_GENES}")
else
    GENES=("${ALL_GENES[@]}")
fi

echo "============================================================" | tee "$LOG_FILE"
echo "BATCH MINIPROT VALIDATION - DUPLICATIONS" | tee -a "$LOG_FILE"
echo "Started: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "Gene file: $GENE_FILE" | tee -a "$LOG_FILE"
echo "Total genes in file: $TOTAL_GENES" | tee -a "$LOG_FILE"
echo "Processing: ${#GENES[@]} genes" | tee -a "$LOG_FILE"
echo "Species: ${#ALL_SPECIES[@]}" | tee -a "$LOG_FILE"
echo "Results: $RESULTS_DIR" | tee -a "$LOG_FILE"
echo "Min identity: $MIN_IDENTITY%" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Track progress
COMPLETED=0
FAILED=0

for GENE in "${GENES[@]}"; do
    ((COMPLETED++)) || true

    echo "[${COMPLETED}/${#GENES[@]}] Processing: $GENE" | tee -a "$LOG_FILE"
    echo "  Time: $(date +%H:%M:%S)" | tee -a "$LOG_FILE"

    OUTPUT_CSV="$RESULTS_DIR/${GENE}_dup_validation.csv"

    if python3 validate_duplications_miniprot.py \
        --genes "$GENE" \
        --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
        --species "$SPECIES_STR" \
        --output "$OUTPUT_CSV" \
        --genome-dir "$GENOME_DIR" \
        --min-n50 5000000 \
        --min-identity "$MIN_IDENTITY" \
        --keep-genomes 2>&1 | tee -a "$LOG_FILE"; then

        if [[ -f "$OUTPUT_CSV" ]]; then
            echo "  ✓ Saved: ${GENE}_dup_validation.csv" | tee -a "$LOG_FILE"
        fi
    else
        echo "  ✗ FAILED: $GENE" | tee -a "$LOG_FILE"
        ((FAILED++)) || true
    fi

    # Storage check every 10 genes
    if (( COMPLETED % 10 == 0 )); then
        AVAIL=$(df -h . | tail -1 | awk '{print $4}')
        echo "  [Storage: $AVAIL available]" | tee -a "$LOG_FILE"
    fi

    echo "" | tee -a "$LOG_FILE"
done

echo "============================================================" | tee -a "$LOG_FILE"
echo "DUPLICATION VALIDATION COMPLETE" | tee -a "$LOG_FILE"
echo "Finished: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "Completed: $COMPLETED genes" | tee -a "$LOG_FILE"
echo "Failed: $FAILED genes" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Results:" | tee -a "$LOG_FILE"
ls -la "$RESULTS_DIR"/*.csv 2>/dev/null | wc -l | xargs -I {} echo "  {} CSV files generated" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
df -h . | tail -1 | tee -a "$LOG_FILE"
