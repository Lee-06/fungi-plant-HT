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
	
    if [ ! -f "${plantgenome}.nsq" ]; then
        echo "Building DB for $plantfilename..."
        makeblastdb -in "$plantgenome" -dbtype nucl
    fi

    for fungigenome in "$FUNGI_DIR"/*.fasta; do
        fungifilename=$(basename "$fungigenome")
        
        echo "Blasting $plantfilename vs $fungifilename..."
        # Note: Using hs-blastn here is highly recommended over blastn for speed
        hs-blastn -db "$plantgenome" -query "$fungigenome" \
            -num_threads 20 -evalue 1e-20 -outfmt 6 \
            -out blastresults/${plantfilename}_VS_${fungifilename}.blast
    done
done
