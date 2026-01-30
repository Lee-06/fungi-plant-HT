#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# ─── ARGUMENT PARSING ───────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Check distribution of candidates across plant genomes.")
parser.add_argument("-s", "--selected_fasta", required=True, help="Directory containing extracted FASTA sequences (from Script 3)")
parser.add_argument("-p", "--plant_genomes", required=True, help="Directory containing all plant genome FASTAs")
parser.add_argument("-j", "--jobs", type=int, default=4, help="Number of parallel genomes to process at once")
parser.add_argument("-t", "--threads", type=int, default=8, help="BLAST threads per job")

args = parser.parse_args()

SELECTED_SEQUENCES_DIR = args.selected_fasta
PLANT_GENOMES_DIR = args.plant_genomes
MAX_PARALLEL_JOBS = args.jobs
THREADS_PER_JOB = args.threads

# Outputs
PLANT_BLAST_RESULTS_DIR = "plant_blast_outputs"
PLANT_COMBINED_RESULTS = "plant_alignment_results.tsv"

# BLAST Parameters
OUTFMT = "6 qseqid sseqid pident length evalue bitscore"
EVALUE_THRESHOLD = 1e-20

# ─── SETUP ──────────────────────────────────────────────────────────────────────
if not os.path.exists(PLANT_BLAST_RESULTS_DIR):
    os.makedirs(PLANT_BLAST_RESULTS_DIR)

# ─── FUNCTIONS ──────────────────────────────────────────────────────────────────
def check_blast_db(fasta_path):
    """Checks if BLAST DB exists for a fasta, creates it if not."""
    # Check for .nsq or .nal or .nin
    if not (os.path.exists(fasta_path + ".nsq") or os.path.exists(fasta_path + ".nin")):
        # Ensure we don't have a race condition if multiple jobs try to make the same DB
        # Ideally, DBs are pre-built by Script 1. We assume they exist or try to build.
        cmd = ["makeblastdb", "-in", fasta_path, "-dbtype", "nucl", "-out", fasta_path]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_blast_job(plant_fasta_path, query_fasta_path, output_path):
    """Runs blastn for a specific plant genome against the query sequences."""
    try:
        # ensure db exists
        check_blast_db(plant_fasta_path)
        
        cmd = [
            "blastn",
            "-db", plant_fasta_path,
            "-query", query_fasta_path,
            "-outfmt", OUTFMT,
            "-evalue", str(EVALUE_THRESHOLD),
            "-num_threads", str(THREADS_PER_JOB),
            "-out", output_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return f"Error blasting {plant_fasta_path}: {result.stderr}"
        return None  # Success
    except Exception as e:
        return str(e)

# ─── MAIN EXECUTION ─────────────────────────────────────────────────────────────
def main():
    # 1. Consolidate all 'selected' sequences into one master query file
    # This is much more efficient than blasting many small fasta files individually.
    master_query_file = "all_candidates_query.fasta"
    print("[INFO] Consolidating candidate sequences into master query...")
    
    with open(master_query_file, 'w') as outfile:
        found_any = False
        for filename in os.listdir(SELECTED_SEQUENCES_DIR):
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                path = os.path.join(SELECTED_SEQUENCES_DIR, filename)
                with open(path, 'r') as infile:
                    outfile.write(infile.read())
                    found_any = True
    
    if not found_any:
        sys.exit(f"[ERROR] No fasta files found in {SELECTED_SEQUENCES_DIR}")

    # 2. Prepare Jobs
    plant_files = [
        os.path.join(PLANT_GENOMES_DIR, f) 
        for f in os.listdir(PLANT_GENOMES_DIR) 
        if f.endswith(('.fasta', '.fa', '.fna'))
    ]
    
    print(f"[INFO] Found {len(plant_files)} plant genomes to check against.")
    
    # 3. Run BLAST in Parallel
    print(f"[INFO] Running BLAST with {MAX_PARALLEL_JOBS} parallel jobs...")
    
    with ProcessPoolExecutor(max_workers=MAX_PARALLEL_JOBS) as executor:
        futures = {}
        for plant_path in plant_files:
            plant_name = os.path.basename(plant_path)
            out_file = os.path.join(PLANT_BLAST_RESULTS_DIR, f"{plant_name}.tsv")
            
            # Skip if already done (resume capability)
            if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
                continue
                
            future = executor.submit(run_blast_job, plant_path, master_query_file, out_file)
            futures[future] = plant_name

        for future in as_completed(futures):
            plant_name = futures[future]
            error = future.result()
            if error:
                print(f"[ERROR] {plant_name}: {error}", file=sys.stderr)
            else:
                # Optional: Print progress dots
                print(".", end="", flush=True)

    print("\n[INFO] BLAST processing complete.")

    # 4. Combine Results
    print("[INFO] Combining results...")
    all_plant_hits = []
    columns_list = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]

    for filename in os.listdir(PLANT_BLAST_RESULTS_DIR):
        if not filename.endswith(".tsv"):
            continue
            
        file_path = os.path.join(PLANT_BLAST_RESULTS_DIR, filename)
        
        # Determine plant name from filename
        plant_genome_name = filename.replace(".tsv", "") # Simplification
        
        try:
            # Check for empty files
            if os.path.getsize(file_path) == 0:
                continue

            df = pd.read_csv(file_path, sep="\t", names=columns_list)
            df["plant_genome"] = plant_genome_name
            all_plant_hits.append(df)
        except Exception as e:
            print(f"[WARNING] Could not read {filename}: {e}")

    if all_plant_hits:
        plant_results_df = pd.concat(all_plant_hits, ignore_index=True)
        
        # Save raw combined results
        plant_results_df.to_csv(PLANT_COMBINED_RESULTS, sep="\t", index=False)
        print(f"[SUCCESS] Combined plant hits saved to {PLANT_COMBINED_RESULTS}")
        print(f"[INFO] Total hits found: {len(plant_results_df)}")
        
        # (The comparison logic happens in Script 5, so we stop here)
    else:
        print("[WARNING] No hits found in any plant genome.")

    # Clean up temp file
    if os.path.exists(master_query_file):
        os.remove(master_query_file)

if __name__ == "__main__":
    main()
