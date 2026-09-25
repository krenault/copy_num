#!/bin/bash
# Run gene validation using NCBI BLAST
#
# Usage:
#   ./run_validation.sh RAB40B "Sciurus carolinensis" 10
#   ./run_validation.sh NKG7 "Myotis davidii" 8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check for required arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 GENE_SYMBOL FOCAL_SPECIES [N_NEIGHBORS]"
    echo ""
    echo "Examples:"
    echo "  $0 RAB40B 'Sciurus carolinensis' 10"
    echo "  $0 NKG7 'Myotis davidii' 8"
    exit 1
fi

GENE="$1"
FOCAL_SPECIES="$2"
N_NEIGHBORS="${3:-10}"

echo "========================================"
echo "Gene Validation via NCBI BLAST"
echo "========================================"
echo "Gene: $GENE"
echo "Focal species: $FOCAL_SPECIES"
echo "N neighbors: $N_NEIGHBORS"
echo "========================================"

cd "$SCRIPT_DIR/../.."

# Check if biopython is installed
python3 -c "import Bio" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Installing biopython..."
    pip3 install biopython requests --quiet
fi

# Run the validation
python3 "$SCRIPT_DIR/validate_gene.py" \
    --gene "$GENE" \
    --focal-species "$FOCAL_SPECIES" \
    --n-neighbors "$N_NEIGHBORS"
