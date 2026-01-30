#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

parser = argparse.ArgumentParser(description="Build Phylogenetic Trees for HGT Validation")

parser.add_argument("-i", "--input", required=True, help="Input FASTA file (e.g., hgt_filtered.fasta)")
parser.add_argument("-db", "--database", required=True, help="Path to local BLAST database (nt or custom combined)")
parser.add_argument("-o", "--outdir", default="phylogenies", help="Output directory for trees")
parser.add_argument("-t", "--threads", default=4, help="Number of threads for MAFFT/IQ-TREE")
parser.add_argument("--max_hits", default=50, type=int, help="Max homologs to retrieve per candidate")

args = parser.parse_args()

def check_tool(name):
    from shutil import which
    if which(name) is None:
        sys.exit(f"Error: {name} is not installed or not in your PATH.")

check_tool("mafft")
check_tool("iqtree")
check_tool("blastn")

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

def run_blast_local(query_seq, db_path, out_file, threads):
    """Runs blastn locally against a specified database."""
    cmd = [
        "blastn", "-db", db_path, "-query", "-", 
        "-outfmt", "6 sseqid sseq", 
        "-max_target_seqs", str(args.max_hits), 
        "-num_threads", str(threads)
    ]
    
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input=query_seq)
    
    if process.returncode != 0:
        print(f"Error in BLAST: {stderr}")
        return []
    
    hits = []
    for line in stdout.strip().split("\n"):
        if line:
            parts = line.split("\t")
            if len(parts) >= 2:
                hits.append((parts[0], parts[1]))
    return hits

def run_mafft(input_fasta, output_aln, threads):
    cmd = ["mafft", "--thread", str(threads), "--auto", input_fasta]
    with open(output_aln, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL)

def run_iqtree(input_aln, threads):
    # -bb 1000 = UltraFast Bootstrap
    cmd = ["iqtree", "-s", input_aln, "-bb", "1000", "-nt", str(threads), "-quiet"]
    subprocess.run(cmd)

# Main Execution
for record in SeqIO.parse(args.input, "fasta"):
    candidate_id = record.id.replace(":", "_").replace("|", "_")
    print(f"--> Processing candidate: {candidate_id}")
    
    # 1. Fetch Homologs
    homologs = run_blast_local(str(record.seq), args.database, None, args.threads)
    
    if not homologs:
        print(f"    No homologs found for {candidate_id}. Skipping tree.")
        continue

    # 2. Create Unaligned Fasta (Candidate + Homologs)
    fasta_path = os.path.join(args.outdir, f"{candidate_id}_homologs.fasta")
    with open(fasta_path, "w") as f:
        # Write Original Candidate
        f.write(f">{candidate_id}_CANDIDATE\n{str(record.seq)}\n")
        # Write Homologs
        seen = set()
        for h_id, h_seq in homologs:
            if h_id not in seen:
                f.write(f">{h_id}\n{h_seq}\n")
                seen.add(h_id)
    
    # 3. Align
    aln_path = os.path.join(args.outdir, f"{candidate_id}.aln")
    run_mafft(fasta_path, aln_path, args.threads)
    
    # 4. Tree
    run_iqtree(aln_path, args.threads)
    print(f"    Tree generated: {aln_path}.treefile")

print("Done. Trees are in", args.outdir)
