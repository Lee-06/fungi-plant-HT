#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
import os 
import pandas as pd

# Path to the directory containing your .blast files
blast_dir = "[TOREPLACE]"  # Replace with your directory path
output_file = "filtered_blast_results_with_fungi.tsv"

# Define filtering thresholds
identity_threshold = 80
alignment_length_threshold = 500

columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

filtered_results = pd.DataFrame()

for filename in os.listdir(blast_dir):
    if filename.endswith(".blast"):
        file_path = os.path.join(blast_dir, filename)
        fungi_genome = filename.replace(".blast", "")
        blast_results = pd.read_csv(file_path, sep="\t", names=columns)
        blast_results["fungi_genome"] = fungi_genome
        filtered = blast_results[(blast_results["pident"] >= identity_threshold) & (blast_results["length"] >= alignment_length_threshold)]
        filtered_results = pd.concat([filtered_results, filtered], ignore_index=True)

filtered_results.to_csv(output_file, sep="\t", index=False)
