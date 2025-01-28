#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adapted script to process multiple FASTA files instead of a single one.
@author: Lee Mariault
"""

import os
import pandas as pd
from Bio import SeqIO

plant_genomes_dir = "[TOREPLACE]"  # Replace with the directory containing the plant genome FASTAs
filtered_results_file = "filtered_blast_results_with_fungi.tsv"
filtered_results = pd.read_csv(filtered_results_file, sep="\t")
selected_sseqids = set(filtered_results["sseqid"])

for genome_file in os.listdir(plant_genomes_dir):
    if genome_file.lower().endswith((".fasta", ".fa", ".fna")):
        fasta_path = os.path.join(plant_genomes_dir, genome_file)
        base_name = os.path.splitext(genome_file)[0]
        output_fasta = f"selected_{base_name}.fasta"
        with open(output_fasta, "w") as output_handle:
            for record in SeqIO.parse(fasta_path, "fasta"):
                if record.id in selected_sseqids:
                    SeqIO.write(record, output_handle, "fasta")
