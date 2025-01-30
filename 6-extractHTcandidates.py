#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
from Bio import SeqIO
import pandas as pd

INPUT_TSV = "hgt_candidates.tsv"
GENOME_OF_INTEREST_FASTA = "[TOREPLACE]" #Fill in with the genome you are interested in (or create a loop)
OUTPUT_FASTA = "hgt_candidates.fasta"

genome_seqs = SeqIO.to_dict(SeqIO.parse(GENOME_OF_INTEREST_FASTA, "fasta"))

hgt_df = pd.read_csv(INPUT_TSV, sep="\t")

with open(OUTPUT_FASTA, "w") as out_fasta:
    for _, row in hgt_df.iterrows():
        chrom = str(row["sseqid_fungi"]).strip()
        start = int(row["sstart_fungi"])
        end = int(row["send_fungi"])

        if chrom in genome_seqs:
            seq_record = genome_seqs[chrom]
            
            if start > end:
                start, end = end, start
            
            extracted_seq = seq_record.seq[start-1:end]
            out_fasta.write(f">{chrom}:{start}-{end}\n{extracted_seq}\n")
        else:
            print(f"[WARNING] Chromosome/Scaffold {chrom} not found in genome FASTA!")
