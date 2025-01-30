#!/bin/bash

cd-hit -i hgt_candidates.fasta -o hgt_clusters.fasta -c 0.9 -n 5 -d 0 -T 8 -M 16000

#eggnog is required for the following step
emapper.py -i hgt_clusters.fasta \
           --output hgt_annotations \
           --cpu 8 \
           --tax_scope 2759 \
           --dmnd_db /media/lee/bigdata/eggnog_database_files/eggnog_proteins.dmnd \
           -m diamond \
           --data_dir /media/lee/bigdata/eggnog_database_files \
           --translate
