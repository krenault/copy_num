#!/bin/bash

# Base URL for TOGA human hg38 reference data
BASE_URL="https://genome.senckenberg.de/download/TOGA/human_hg38_reference"

# Output directory
OUTPUT_DIR="file.path(ROOT, "data")"
mkdir -p "$OUTPUT_DIR"

# File to store unique species names
SPECIES_LIST_FILE="unique_species_list.txt"
> "$SPECIES_LIST_FILE"  # Create or truncate the file

# Function to safely remove directory if empty (macOS compatible)
remove_if_empty() {
    if [ -d "$1" ] && [ -z "$(ls -A "$1")" ]; then
        rmdir "$1"
        echo "Removed empty directory: $1"
    fi
}

# Get list of all order directories (first level)
echo "Discovering order directories..."
ORDER_DIRS=$(curl -sk "$BASE_URL/" | 
             grep -o 'href="[^"]*/"' | 
             grep -v 'Parent\|Multiple\|chain\|overview\|Matrix\|hg38' | 
             sed 's/href="//g; s/"//g; s/\///g')

if [ -z "$ORDER_DIRS" ]; then
    echo "Error: No order directories found"
    exit 1
fi

# Process each order directory
for order in $ORDER_DIRS; do
    echo "Processing order: $order"
    mkdir -p "$OUTPUT_DIR/$order"

    # Get list of species directories (second level)
    echo "Getting species for $order..."
    SPECIES_DIRS=$(curl -sk "$BASE_URL/$order/" | 
                   grep -o 'href="[^"]*/"' | 
                   grep -v 'Parent' | 
                   sed 's/href="//g; s/"//g; s/\///g')

    for species in $SPECIES_DIRS; do
        echo "Processing species: $species"
        
        # Extract species name (remove any IDs or suffixes)
        SPECIES_NAME=$(echo "$species" | sed 's/__.*$//' | tr '_' ' ')
        
        # Add species name to the list if not already present
        if ! grep -q "^$SPECIES_NAME$" "$SPECIES_LIST_FILE"; then
            echo "$SPECIES_NAME" >> "$SPECIES_LIST_FILE"
            echo "Added $SPECIES_NAME to species list"
        fi
        
        mkdir -p "$OUTPUT_DIR/$order/$species"
        
        # Get all files in the species directory
        SPECIES_FILES=$(curl -sk "$BASE_URL/$order/$species/" | 
                        grep -o 'href="[^"]*"' | 
                        sed 's/href="//g; s/"//g' | 
                        grep -v 'Parent\|Directory')
        
        # Filter for files ending with loss_summ_data.tsv
        LOSS_FILES=$(echo "$SPECIES_FILES" | grep ".*orthologsClassification\.tsv.gz$")
        
        FOUND_FILES=false
        
        if [ -n "$LOSS_FILES" ]; then
            for file in $LOSS_FILES; do
                LOSS_FILE_URL="${BASE_URL}/${order}/${species}/${file}"
                OUTPUT_FILE="${OUTPUT_DIR}/${order}/${species}/${file}"
                
                echo "Found matching file: $file"
                
                # Check if file already exists and is not empty
                if [ -s "$OUTPUT_FILE" ]; then
                    echo "File already exists and is not empty, skipping: $OUTPUT_FILE"
                    FOUND_FILES=true
                    continue
                fi
                
                echo "Downloading from: $LOSS_FILE_URL"
                if curl -sk -f "$LOSS_FILE_URL" > "$OUTPUT_FILE" 2>/dev/null; then
                    if [ -s "$OUTPUT_FILE" ]; then
                        echo "Successfully saved to: $OUTPUT_FILE"
                        echo "First line: $(head -1 "$OUTPUT_FILE" | cut -c1-80)..."
                        FOUND_FILES=true
                    else
                        echo "Warning: Downloaded file is empty"
                        rm -f "$OUTPUT_FILE"
                    fi
                else
                    echo "Failed to download: $file"
                    rm -f "$OUTPUT_FILE"
                fi
            done
        else
            echo "No files ending with loss_summ_data.tsv found for $species"
        fi
        
        if [ "$FOUND_FILES" = false ]; then
            # Remove directory if empty
            remove_if_empty "$OUTPUT_DIR/$order/$species"
        fi
        
        echo ""
    done
    
    # Remove order directory if empty
    remove_if_empty "$OUTPUT_DIR/$order"
done

# Sort the species list alphabetically
sort -o "$SPECIES_LIST_FILE" "$SPECIES_LIST_FILE"

# Count the number of unique species
SPECIES_COUNT=$(wc -l < "$SPECIES_LIST_FILE")
echo "All *loss_summ_data.tsv files downloaded to: $OUTPUT_DIR"
echo "Created list of $SPECIES_COUNT unique species in $SPECIES_LIST_FILE"
