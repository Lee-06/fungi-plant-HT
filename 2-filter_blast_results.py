#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
import os 
import pandas as pd

# Path to the directory containing your .blast files
blast_dir = "[TOREPLACE]"  # Replace with your directory path
# Directories containing the .fasta.fai index files for fungi and plants
fungi_fai_dir = "[TOREPLACE]"  # Replace with fungi FASTA index directory
plant_fai_dir = "[TOREPLACE]"  # Replace with plant FASTA index directory
output_file = "filtered_blast_results_with_fungi.tsv"

# Define filtering thresholds
identity_threshold = 80
alignment_length_threshold = 500

scaffold_length_threshold = 20000  # 20 kb

columns = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def load_fai(fai_path):
    lengths = {}
    if os.path.exists(fai_path):
        with open(fai_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
    return lengths

filtered_results = pd.DataFrame()

for filename in os.listdir(blast_dir):
    if filename.endswith(".blast"):
        file_path = os.path.join(blast_dir, filename)
        base = filename[:-6]  # remove .blast
        if "_VS_" in base:
            plant_name, fungi_name = base.split("_VS_", 1)
        else:
            # fallback if naming convention is different
            fungi_name = base
            plant_name = ""

        fungi_fai = os.path.join(fungi_fai_dir, fungi_name + ".fai")
        plant_fai = os.path.join(plant_fai_dir, plant_name + ".fai") if plant_name else None

        fungi_lengths = load_fai(fungi_fai)
        plant_lengths = load_fai(plant_fai) if plant_fai else {}
        blast_results = pd.read_csv(file_path, sep="\t", names=columns)
        blast_results["fungi_genome"] = fungi_name
        filtered = blast_results[
            (blast_results["pident"] >= identity_threshold)
            & (blast_results["length"] >= alignment_length_threshold)
            & (blast_results["qseqid"].map(lambda x: fungi_lengths.get(x, 0)) >= scaffold_length_threshold)
            & (blast_results["sseqid"].map(lambda x: plant_lengths.get(x, 0)) >= scaffold_length_threshold)
        ]
        filtered_results = pd.concat([filtered_results, filtered], ignore_index=True)

filtered_results.to_csv(output_file, sep="\t", index=False)
