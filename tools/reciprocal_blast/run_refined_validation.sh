#!/bin/bash
# Validate refined longevity gene candidates with miniprot
# All candidates are known longevity genes with proper longevity associations

set -e
cd "$(dirname "$0")"

GENOME_DIR="genomes"
RESULTS_DIR="results"
MIN_IDENTITY=80

mkdir -p "$GENOME_DIR" "$RESULTS_DIR"

echo "============================================================"
echo "REFINED LONGEVITY CANDIDATE VALIDATION WITH MINIPROT"
echo "============================================================"
echo ""
echo "All candidates are:"
echo "  - Known longevity genes from de Magalhães 2023"
echo "  - Either rare clustered losses in long-lived species"
echo "  - Or convergent changes across diverged long-lived species"
echo ""

# AFFECTED SPECIES (where losses/duplications occurred)
# Plus relatives and outgroups for comparison

# Key species for validation
VALIDATION_SPECIES=(
    # Affected Chiroptera (bats) - many candidates
    "Myotis_brandtii"        # CDC42 loss (MLres=0.764!)
    "Antrozous_pallidus"     # PDGFRB loss
    "Myotis_septentrionalis" # GSTP1 loss
    "Eidolon_helvum"         # APOE loss

    # Bat controls (no loss expected)
    "Myotis_myotis"
    "Pteropus_vampyrus"

    # Affected Primates
    "Gorilla_gorilla"        # FOXO1 duplication
    "Cheirogaleus_medius"    # TXN loss
    "Eulemur_fulvus"         # CHEK2 loss
    "Rhinopithecus_bieti"    # ZMPSTE24 loss

    # Primate controls
    "Homo_sapiens"
    "Pan_troglodytes"
    "Macaca_mulatta"

    # Affected Rodentia
    "Sciurus_carolinensis"   # ZMPSTE24 loss, HTRA2 loss

    # Rodent controls
    "Mus_musculus"
    "Rattus_norvegicus"

    # Affected Cingulata
    "Tolypeutes_matacus"     # PRKCD loss (armadillo - like CGAS!)

    # Carnivore controls
    "Panthera_tigris"
    "Canis_lupus_familiaris"

    # Other long-lived outgroups
    "Elephas_maximus"
    "Equus_caballus"
)

SPECIES_STR=$(IFS=,; echo "${VALIDATION_SPECIES[*]}")

echo "Species for validation (${#VALIDATION_SPECIES[@]} total):"
echo "  Affected + relatives + outgroups"
echo ""

# TOP REFINED CANDIDATES

echo "============================================================"
echo "PRIORITY 2: RARE CLUSTERED LOSSES IN LONG-LIVED (like CGAS)"
echo "============================================================"

# Single losses in definitively long-lived species
RARE_LOSS_GENES="CDC42,PDGFRB,TXN,PRKCD,CHEK2,PDGFRA"

echo ""
echo "Genes: $RARE_LOSS_GENES"
echo "Pattern: Single loss in species with high MLres (0.2-0.8)"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "$RARE_LOSS_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_rare_loss_longevity.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 2: RARE CLUSTERED DUPLICATIONS"
echo "============================================================"

RARE_DUP_GENES="FOXO1"

echo ""
echo "Genes: $RARE_DUP_GENES"
echo "FOXO1 duplication in Gorilla (FOXO transcription factor!)"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "$RARE_DUP_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_rare_dup_longevity.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 3: CONVERGENT LOSSES IN LONG-LIVED"
echo "============================================================"

# Multiple losses across clades, all in long-lived species
CONVERGENT_LOSS_GENES="ZMPSTE24,GSTP1,HTRA2,APOE,CDK7,TRPV1"

echo ""
echo "Genes: $CONVERGENT_LOSS_GENES"
echo "Pattern: Losses in 2-10 species across multiple orders"
echo "         >50% of affected species are long-lived (MLres > 0)"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "$CONVERGENT_LOSS_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_loss_longevity.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "PRIORITY 3: CONVERGENT DUPLICATIONS IN LONG-LIVED"
echo "============================================================"

CONVERGENT_DUP_GENES="APEX1,GCLM,CDKN1A,MTOR"

echo ""
echo "Genes: $CONVERGENT_DUP_GENES"
echo "Pattern: Duplications in long-lived species across clades"
echo ""

python3 validate_candidates_miniprot.py \
    --genes "$CONVERGENT_DUP_GENES" \
    --toga-matrix "../../data/duplication_dollo/All_Species_Gene_Duplication_Binary.tsv" \
    --species "$SPECIES_STR" \
    --output "$RESULTS_DIR/miniprot_convergent_dup_longevity.csv" \
    --genome-dir "$GENOME_DIR" \
    --min-n50 5000000 \
    --min-identity "$MIN_IDENTITY" \
    --keep-genomes 2>&1 || echo "Some genes failed - continuing"

echo ""
echo "============================================================"
echo "VALIDATION COMPLETE"
echo "============================================================"
echo ""
echo "Results saved to:"
ls -la "$RESULTS_DIR"/miniprot_*_longevity.csv 2>/dev/null || echo "  (checking...)"
echo ""
echo "Key interpretations:"
echo "  - If miniprot shows 0 loci where TOGA says lost -> TRUE LOSS (validate)"
echo "  - If miniprot shows 1 locus where TOGA says lost -> FALSE POSITIVE"
echo "  - If miniprot shows >1 loci where TOGA says dup -> TRUE DUPLICATION"
echo ""
echo "Top candidates to highlight in presentation:"
echo "  1. CDC42 - loss in long-lived bat (MLres=0.76) - cell cycle"
echo "  2. ZMPSTE24 - losses in snub-nosed monkey + squirrel - PROGERIA"
echo "  3. APOE - losses across primates/bats - ALZHEIMER'S"
echo "  4. FOXO1 - duplication in Gorilla - FOXO longevity pathway"
