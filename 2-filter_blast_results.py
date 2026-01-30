#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
import os
import sys
import argparse
import pandas as pd

# ─── ARGUMENT PARSING ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Filter BLAST results for HT candidates.")
parser.add_argument("--blast_dir", required=True, help="Directory containing .blast output files")
parser.add_argument("--fungi_fai", required=True, help="Directory containing Fungi .fasta.fai index files")
parser.add_argument("--plant_fai", required=True, help="Directory containing Plant .fasta.fai index files")
parser.add_argument("--output", default="filtered_blast_results_with_fungi.tsv", help="Output TSV filename")

args = parser.parse_args()

BLAST_DIR = args.blast_dir
FUNGI_FAI_DIR = args.fungi_fai
PLANT_FAI_DIR = args.plant_fai
OUTPUT_FILE = args.output

# ─── CONFIGURATION ──────────────────────────────────────────────────────────────
# Thresholds defined in the thesis context
IDENTITY_THRESHOLD = 80
ALIGNMENT_LENGTH_THRESHOLD = 500
SCAFFOLD_LENGTH_THRESHOLD = 20000  # 20 kb (filters out short/fragmented scaffolds)

COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]

# ─── HELPER FUNCTIONS ───────────────────────────────────────────────────────────
def load_fai(fai_path):
    """Parses a .fai file and returns a dictionary of {seq_id: length}."""
    lengths = {}
    if os.path.exists(fai_path):
        try:
            with open(fai_path, 'r') as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        lengths[parts[0]] = int(parts[1])
        except Exception as e:
            print(f"[WARNING] Could not read FAI file {fai_path}: {e}", file=sys.stderr)
    return lengths

# ─── MAIN EXECUTION ─────────────────────────────────────────────────────────────
def main():
    if not os.path.isdir(BLAST_DIR):
        sys.exit(f"[ERROR] BLAST directory not found: {BLAST_DIR}")

    print(f"[INFO] Starting filtering process...")
    print(f"[INFO] Thresholds: Identity>={IDENTITY_THRESHOLD}%, Len>={ALIGNMENT_LENGTH_THRESHOLD}bp, Scaffold>={SCAFFOLD_LENGTH_THRESHOLD}bp")

    filtered_frames = []  # List to store valid dataframes (much faster than append)
    file_count = 0

    for filename in os.listdir(BLAST_DIR):
        if not filename.endswith(".blast"):
            continue
            
        file_path = os.path.join(BLAST_DIR, filename)
        base = filename[:-6]  # remove .blast suffix

        # Parse filename to get genome names
        if "_VS_" in base:
            plant_name, fungi_name = base.split("_VS_", 1)
        else:
            # Fallback logic
            fungi_name = base
            plant_name = ""

        # Paths to index files
        fungi_fai_path = os.path.join(FUNGI_FAI_DIR, fungi_name + ".fasta.fai")
        # Handle inconsistent extensions if necessary (.fai or .fasta.fai)
        if not os.path.exists(fungi_fai_path):
             fungi_fai_path = os.path.join(FUNGI_FAI_DIR, fungi_name + ".fai")

        # Load scaffold lengths
        fungi_lengths = load_fai(fungi_fai_path)
        
        # (Optional: Load plant lengths if you need to filter by plant scaffold size too)
        # plant_fai_path = os.path.join(PLANT_FAI_DIR, plant_name + ".fai")
        
        if not fungi_lengths:
            # If we can't verify scaffold length, we might skip or warn. 
            # Here we warn but proceed (filtering will fail for missing IDs).
            # print(f"[DEBUG] No FAI found for {fungi_name}", file=sys.stderr)
            pass

        try:
            # Read BLAST result
            # Check if file is empty first to avoid pandas errors
            if os.path.getsize(file_path) == 0:
                continue

            blast_results = pd.read_csv(file_path, sep="\t", names=COLUMNS)
            blast_results["fungi_genome"] = fungi_name
            
            # Apply Filters
            # 1. Identity & Alignment Length
            mask = (blast_results["pident"] >= IDENTITY_THRESHOLD) & \
                   (blast_results["length"] >= ALIGNMENT_LENGTH_THRESHOLD)
            
            temp_filtered = blast_results[mask]

            if temp_filtered.empty:
                continue

            # 2. Scaffold Length (Map qseqid to length dictionary)
            # This ensures the fungal hit is on a substantial scaffold, not a tiny contig
            temp_filtered = temp_filtered[
                temp_filtered["qseqid"].map(lambda x: fungi_lengths.get(x, 0)) >= SCAFFOLD_LENGTH_THRESHOLD
            ]

            if not temp_filtered.empty:
                filtered_frames.append(temp_filtered)
        
        except Exception as e:
            print(f"[ERROR] Processing {filename}: {e}", file=sys.stderr)
        
        file_count += 1
        if file_count % 1000 == 0:
            print(f"[INFO] Processed {file_count} files...")

    # Concatenate all results
    if filtered_frames:
        final_df = pd.concat(filtered_frames, ignore_index=True)
        final_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
        print(f"[SUCCESS] Filtered results saved to {OUTPUT_FILE}")
        print(f"[INFO] Total hits kept: {len(final_df)}")
    else:
        print("[WARNING] No hits passed the filters. Output file not created.")

if __name__ == "__main__":
    main()
