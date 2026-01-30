#!/bin/bash

if [[ -z "$1" ]]; then
    echo "Error: Missing EggNOG database path."
    echo "Usage: $0 <path_to_eggnog_data_dir>"
    exit 1
fi

EGGNOG_DATA_DIR="$1"

# Run CD-HIT
cd-hit -i hgt_candidates.fasta -o hgt_clusters.fasta -c 0.9 -n 5 -d 0 -T 8 -M 16000

# Run EggNOG Mapper using the provided path
emapper.py -i hgt_clusters.fasta \
           --output hgt_annotations \
           --cpu 8 \
           --tax_scope 2759 \
           --dmnd_db "${EGGNOG_DATA_DIR}/eggnog_proteins.dmnd" \
           -m diamond \
           --data_dir "${EGGNOG_DATA_DIR}" \
           --translate
