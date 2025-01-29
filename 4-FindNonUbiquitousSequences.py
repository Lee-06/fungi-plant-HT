#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""

import os
import subprocess
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

# Required inputs
SELECTED_SEQUENCES_FASTA = "[TOREPLACE]"   # Directory with extracted fasta sequences
PLANT_GENOMES_DIR        = "[TOREPLACE]"    # Directory with plant genome FASTAs
FUNGI_RESULTS_FILE       = "filtered_blast_results_with_fungi.tsv"

# Intermediate outputs
PLANT_BLAST_RESULTS_DIR = "plant_blast_outputs"         # Directory to store each plant BLAST result
PLANT_COMBINED_RESULTS  = "plant_alignment_results.tsv"  # All plant BLAST results combined
FUNGI_VS_PLANT_FILE     = "fungi_vs_plant_comparison.tsv" # Final comparison file

# BLAST parameters
OUTFMT            = "6 qseqid sseqid pident length evalue bitscore"
EVALUE_THRESHOLD  = 1e-20
MAX_TARGET_SEQS   = 5
#Replace the following according to your computer/cluster specifications
NUM_PARALLEL_JOBS = [TOREPLACE]      # Number of plant genomes to BLAST in parallel
THREADS_PER_JOB   = [TOREPLACE]      # Threads used by each blastn process 

if not os.path.exists(PLANT_BLAST_RESULTS_DIR):
    os.makedirs(PLANT_BLAST_RESULTS_DIR)

# STEP 0: CREATE OR CHECK BLAST DATABASES (SERIAL)

plant_files = [
    f for f in os.listdir(PLANT_GENOMES_DIR)
    if f.lower().endswith((".fasta", ".fa", ".fna"))
]

db_paths = []
for plant_file in plant_files:
    plant_path = os.path.join(PLANT_GENOMES_DIR, plant_file)
    db_name = plant_path
    db_base = os.path.basename(db_name)
    if not any(os.path.exists(db_name + ext) for ext in [".nsq", ".nin", ".ndb"]):
        print(f"[INFO] Creating BLAST DB for {db_base}...")
        makeblastdb_cmd = [
            "makeblastdb",
            "-in", plant_path,
            "-dbtype", "nucl",
            "-out", db_name
        ]
        subprocess.run(makeblastdb_cmd, check=True)
    else:
        print(f"[INFO] BLAST DB for {db_base} already exists. Skipping creation.")

    db_paths.append((plant_path, db_name, db_base))

# STEP 1: GATHER ALL FILTERED (QUERY) FASTA FILES
selected_files = [
    f for f in os.listdir(SELECTED_SEQUENCES_FASTA)
    if f.lower().endswith((".fasta", ".fa", ".fna"))
]

if not selected_files:
    print("[ERROR] No suspect (query) FASTA files found in SUSPECT_SEQUENCES_DIR.")
    exit(1)
else:
    print(f"[INFO] Found {len(selected_files)} suspect FASTA files in {SELECTED_SEQUENCES_FASTA}.")


# STEP 2: PARALLEL BLAST AGAINST EACH PLANT GENOME
def run_blast(query_path, db_info):
    plant_path, db_name, db_base = db_info
    query_base = os.path.splitext(os.path.basename(query_path))[0]
    result_file = os.path.join(PLANT_BLAST_RESULTS_DIR, f"{query_base}_vs_{db_base}.tsv")

    if not os.path.exists(result_file):
        blastn_cmd = [
            "blastn",
            "-query", query_path,
            "-db", db_name,
            "-out", result_file,
            "-outfmt", OUTFMT,
            "-evalue", str(EVALUE_THRESHOLD),
            "-max_target_seqs", str(MAX_TARGET_SEQS),
            "-num_threads", str(THREADS_PER_JOB)
        ]
        subprocess.run(blastn_cmd, check=True)
    else:
        print(f"[INFO] BLAST results for {query_base} vs {db_base} already exist, skipping.")

    return result_file

blast_tasks = []
for query_file in suspect_files:
    query_path = os.path.join(SELECTED_SEQUENCES_FASTA, query_file)
    for db_info in db_paths:
        blast_tasks.append((query_path, db_info))

results_collected = []
with ProcessPoolExecutor(max_workers=NUM_PARALLEL_JOBS) as executor:
    future_to_task = {
        executor.submit(run_blast, query_path, db_info): (query_path, db_info)
        for (query_path, db_info) in blast_tasks
    }
    for future in as_completed(future_to_task):
        query_path, db_info = future_to_task[future]
        query_base = os.path.splitext(os.path.basename(query_path))[0]
        db_base = db_info[2]
        try:
            result_file = future.result()
            results_collected.append(result_file)
            print(f"[INFO] Completed BLAST for {query_base} vs {db_base}; results in {result_file}")
        except Exception as exc:
            print(f"[ERROR] BLAST failed for {query_base} vs {db_base}: {exc}")

# STEP 3: COMBINE AND FIND BEST HITS PER QUERY ACROSS ALL PLANT RESULTS
columns = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
all_plant_hits = []

for res_file in os.listdir(PLANT_BLAST_RESULTS_DIR):
    if not res_file.endswith(".tsv"):
        continue
    file_path = os.path.join(PLANT_BLAST_RESULTS_DIR, res_file)
    plant_genome_name = os.path.splitext(res_file)[0]

    df = pd.read_csv(file_path, sep="\t", names=columns)
    if df.empty:
        continue

    df["result_file"] = res_file
    all_plant_hits.append(df)

if all_plant_hits:
    plant_results_df = pd.concat(all_plant_hits, ignore_index=True)
    best_plant_hits = (
        plant_results_df
        .sort_values("bitscore", ascending=False)
        .groupby("qseqid", as_index=False)
        .first()
    )
    # Save all combined hits to a TSV
    plant_results_df.to_csv(PLANT_COMBINED_RESULTS, sep="\t", index=False)
    print(f"[INFO] Combined plant hits: {plant_results_df.shape[0]} total rows.")
else:
    best_plant_hits = pd.DataFrame(columns=columns + ["plant_genome"])
    print("[WARNING] No hits found in any plant BLAST results.")

best_plant_hits = best_plant_hits.rename(columns={
    "sseqid": "sseqid_plant",
    "pident": "pident_plant",
    "length": "length_plant",
    "evalue": "evalue_plant",
    "bitscore": "bitscore_plant",
    "plant_genome": "best_plant_genome"
})
