#!/bin/bash
# Validate top longevity gene candidates with miniprot
# Based on CGAS-like patterns (rare clustered) and convergent changes in long-lived

set -e
cd "$(dirname "$0")"

GENOME_DIR="genomes"
RESULTS_DIR="results"
MIN_IDENTITY=80

mkdir -p "$GENOME_DIR" "$RESULTS_DIR"

echo "============================================================"
echo "TOP LONGEVITY CANDIDATE VALIDATION WITH MINIPROT"
echo "============================================================"

# Species to validate across (representing major clades)
SPECIES=(
    # Bats (Chiroptera) - key clade for many candidates
    "Myotis_brandtii"
    "Myotis_myotis"
    "Antrozous_pallidus"
    "Eonycteris_spelaea"
    "Pteropus_vampyrus"

    # Long-lived primates
    "Homo_sapiens"
    "Pan_troglodytes"
    "Gorilla_gorilla"
    "Macaca_mulatta"

    # Rodents (short-lived reference)
    "Mus_musculus"
    "Rattus_norvegicus"
    "Sciurus_carolinensis"  # Squirrel - long-lived rodent

    # Carnivores
    "Panthera_tigris"
    "Canis_lupus_familiaris"

    # Other long-lived
    "Elephas_maximus"
    "Balaenoptera_physalus"  # Whale
)

SPECIES_STR=$(IFS=,; echo "${SPECIES[*]}")

# TOP CANDIDATES BY PATTERN

echo ""
echo "============================================================"
echo "PATTERN 1: RARE CLUSTERED LOSSES (like CGAS)"
echo "============================================================"

# These have 1-3 losses in specific clades
RARE_LOSS_GENES="CDC42,PDGFRB,LEPR,CNR1,SPRTN,MLH1,SIRT7"

echo "Genes: $RARE_LOSS_GENES"
python3 validate_candidates_miniprot.py \
    --genes "$RARE_LOSS_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_rare_clustered_loss.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes || echo "Some genes failed"

echo ""
echo "============================================================"
echo "PATTERN 2: CONVERGENT LOSSES IN LONG-LIVED"
echo "============================================================"

# Losses in species from different orders, all long-lived
CONVERGENT_LOSS_GENES="ZMPSTE24,SIRT3,APOE,ATM,GSTP1,HTRA2,TRAP1"

echo "Genes: $CONVERGENT_LOSS_GENES"
python3 validate_candidates_miniprot.py \
    --genes "$CONVERGENT_LOSS_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_loss.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes || echo "Some genes failed"

echo ""
echo "============================================================"
echo "PATTERN 3: CONVERGENT DUPLICATIONS IN LONG-LIVED"
echo "============================================================"

# Duplications in long-lived species across clades
CONVERGENT_DUP_GENES="SIRT6,MTOR,NFE2L2,FOXO1,APEX1,CDKN1A"

echo "Genes: $CONVERGENT_DUP_GENES"
python3 validate_candidates_miniprot.py \
    --genes "$CONVERGENT_DUP_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_dup.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes || echo "Some genes failed"

echo ""
echo "============================================================"
echo "VALIDATION COMPLETE"
echo "============================================================"
echo ""
echo "Results:"
ls -la "$RESULTS_DIR"/miniprot_*.csv 2>/dev/null
echo ""
echo "Summary: Compare miniprot loci to TOGA predictions"
echo "- If miniprot shows 0 loci where TOGA says lost -> TRUE LOSS"
echo "- If miniprot shows 1 locus where TOGA says lost -> FALSE POSITIVE"
echo "- If miniprot shows >1 loci where TOGA says dup -> TRUE DUPLICATION"
