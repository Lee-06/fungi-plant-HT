#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
import pandas as pd

parser = argparse.ArgumentParser(description="Extract FASTA sequences for identified Plant hits.")
parser.add_argument("-i", "--input_tsv", default="filtered_blast_results_with_fungi.tsv", help="Input filtered BLAST results")
parser.add_argument("-p", "--plant_genomes", required=True, help="Directory containing Plant Genome FASTA files")
parser.add_argument("-o", "--outdir", default="selected_sequences", help="Output directory for extracted sequences")

args = parser.parse_args()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

print(f"[INFO] Reading {args.input_tsv}...")
try:
    filtered_results = pd.read_csv(args.input_tsv, sep="\t")
    selected_sseqids = set(filtered_results["sseqid"].astype(str))
except Exception as e:
    exit(f"[ERROR] Could not read input TSV: {e}")

print(f"[INFO] Looking for {len(selected_sseqids)} unique sequences in {args.plant_genomes}...")

found_count = 0
for genome_file in os.listdir(args.plant_genomes):
    if genome_file.lower().endswith((".fasta", ".fa", ".fna")):
        fasta_path = os.path.join(args.plant_genomes, genome_file)
        base_name = os.path.splitext(genome_file)[0]
        output_fasta = os.path.join(args.outdir, f"selected_{base_name}.fasta")
        hits_in_genome = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in selected_sseqids:
                hits_in_genome.append(record)
                selected_sseqids.remove(record.id)
                found_count += 1
        
        if hits_in_genome:
            with open(output_fasta, "w") as output_handle:
                SeqIO.write(hits_in_genome, output_handle, "fasta")

print(f"[SUCCESS] Extracted {found_count} sequences into '{args.outdir}/'.")
if len(selected_sseqids) > 0:
    print(f"[WARNING] {len(selected_sseqids)} sequences were not found in the provided genome directory.")
