#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description="Extract full genomic sequences of HT candidates from Fungal Genomes.")
parser.add_argument("-i", "--input_candidates", required=True, help="Candidate TSV file (e.g., ht_candidates.tsv)")
parser.add_argument("-g", "--fungi_genomes", required=True, help="Directory containing Fungi Genome FASTA files")
parser.add_argument("-o", "--output", default="ht_candidates.fasta", help="Output Multi-FASTA file")

args = parser.parse_args()

df = pd.read_csv(args.input_candidates, sep="\t")
# We need 'fungi_genome' column (added in Script 2) and coordinates
required_cols = ["sseqid_fungi", "sstart_fungi", "send_fungi", "fungi_genome"]
if not all(col in df.columns for col in required_cols):
    exit(f"[ERROR] Input TSV missing one of required columns: {required_cols}")

print(f"[INFO] Extracting {len(df)} sequences...")

grouped = df.groupby("fungi_genome")

with open(args.output, "w") as out_f:
    for genome_name, group in grouped:
        genome_path = os.path.join(args.fungi_genomes, str(genome_name) + ".fasta") 
        if not os.path.exists(genome_path):
             genome_path = os.path.join(args.fungi_genomes, str(genome_name))
        
        if not os.path.exists(genome_path):
            print(f"[WARNING] Genome file for {genome_name} not found in {args.fungi_genomes}. Skipping {len(group)} candidates.")
            continue
            
        print(f"  -> Processing {genome_name}...")
        
        seq_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        
        for _, row in group.iterrows():
            scaffold = str(row["sseqid_fungi"])
            start = int(row["sstart_fungi"])
            end = int(row["send_fungi"])
            
            if start > end: start, end = end, start
            
            if scaffold in seq_dict:
                seq_record = seq_dict[scaffold]
                fragment = seq_record.seq[start-1:end]
                header = f">{row['qseqid']}|{genome_name}|{scaffold}:{start}-{end}"
                out_f.write(f"{header}\n{fragment}\n")
            else:
                print(f"[WARNING] Scaffold {scaffold} not found in {genome_name}")

print(f"[SUCCESS] Extraction complete. Saved to {args.output}")
