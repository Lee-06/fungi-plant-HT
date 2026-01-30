#!/bin/bash

if [[ -z "$1" ]] || [[ -z "$2" ]]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <path_to_plant_genomes> <path_to_fungi_genomes>"
    echo "Example: ./1-BlastWholeGenomes.sh ./genomes/plants ./genomes/fungi"
    exit 1
fi

PLANT_DIR="$1"
FUNGI_DIR="$2"

mkdir -p blastresults

for plantgenome in "$PLANT_DIR"/*.fasta; do
    plantfilename=$(basename "$plantgenome")
    if [ ! -f "${plantgenome}.sys" ]; then
        echo "Building hs-blastn index for $plantfilename..."
        hs-blastn index -i "$plantgenome" -d "$plantgenome"
    fi
    for fungigenome in "$FUNGI_DIR"/*.fasta; do
        fungifilename=$(basename "$fungigenome")
        echo "Blasting $plantfilename vs $fungifilename..."
        hs-blastn align \
            -d "$plantgenome" \
            -q "$fungigenome" \
            -p 20 \
            -f 6 \
            -o "blastresults/${plantfilename}_VS_${fungifilename}.blast"
    done
done
