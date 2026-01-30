#!/usr/bin/env python3
import argparse
import pandas as pd
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter out housekeeping genes based on EggNOG annotations.")
parser.add_argument("--annotations", required=True, help="EggNOG annotation file (.annotations)")
parser.add_argument("--fasta_in", required=True, help="Clustered FASTA input (ht_clusters.fasta)")
parser.add_argument("--fasta_out", default="ht_filtered.fasta", help="Final Filtered FASTA output")

args = parser.parse_args()

# Keywords to exclude
HOUSEKEEPING_KEYWORDS = [
    "ribosomal", "18s", "28s", "5s", "rrna", "rdna", "ribonucleoprotein",
    "translation elongation factor", "mitochondrion", "mitochondrial",
    "cytochrome", "cox1", "nad", "atp6", "atp9", "chloroplast", "plastid",
    "glycolysis", "atp synthase", "nadh dehydrogenase", "oxidoreductase",
    "succinate dehydrogenase", "malate dehydrogenase", "dna polymerase",
    "rna polymerase", "helicase", "topoisomerase", "exonuclease", "ligase",
    "primase", "actin", "tubulin", "kinesin", "dynein", "myosin", "chaperone",
    "heat shock protein", "ubiquitin", "kinase", "phosphatase"
]
keywords_lower = [kw.lower() for kw in HOUSEKEEPING_KEYWORDS]

print("[INFO] Parsing annotations...")
# EggNOG annotation files often have variable headers. We skip comments ##
try:
    df = pd.read_csv(args.annotations, sep="\t", comment="#", header=None)
    # Usually column 0 is Query, column 7 or similar is Description/Free Text.
    # We'll just concat all text columns to be safe for keywords search
except Exception as e:
    exit(f"[ERROR] Could not read annotations: {e}")

# Identify Housekeeping IDs
hk_ids = set()
for index, row in df.iterrows():
    # Convert entire row to string to search for keywords
    row_text = " ".join(row.astype(str)).lower()
    if any(kw in row_text for kw in keywords_lower):
        # Assuming column 0 is the Query ID (standard for EggNOG)
        hk_ids.add(str(row[0]))

print(f"[INFO] Found {len(hk_ids)} sequences matching housekeeping keywords.")

# Write Filtered FASTA
kept_count = 0
excluded_count = 0

with open(args.fasta_out, "w") as out_handle:
    for record in SeqIO.parse(args.fasta_in, "fasta"):
        if record.id not in hk_ids:
            SeqIO.write(record, out_handle, "fasta")
            kept_count += 1
        else:
            excluded_count += 1

print(f"[SUCCESS] Filtered FASTA saved to {args.fasta_out}")
print(f"         Kept: {kept_count} | Removed: {excluded_count}")
