#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
import pandas as pd
from Bio import SeqIO

# Input files
ANNOTATION_FILE = "hgt_annotations.emapper.annotations"
HOUSEKEEPING_BLAST_FILE = "hgt_annotations.emapper.seed_orthologs"
INPUT_FASTA = "hgt_clusters.fasta"
OUTPUT_FASTA = "hgt_filtered.fasta"

# Housekeeping gene keywords to filter
HOUSEKEEPING_KEYWORDS = [
    "ribosomal", "18S", "28S", "5S", "rRNA", "rDNA", "ribonucleoprotein", "translation elongation factor",
    "mitochondrion", "mitochondrial", "cytochrome", "cox1", "nad", "atp6", "atp9", "chloroplast", "plastid",
    "glycolysis", "ATP synthase", "NADH dehydrogenase", "oxidoreductase", "succinate dehydrogenase", "malate dehydrogenase",
    "DNA polymerase", "RNA polymerase", "helicase", "topoisomerase", "exonuclease", "ligase", "primase",
    "actin", "tubulin", "kinesin", "dynein", "myosin",
    "chaperone", "heat shock protein", "ubiquitin", "kinase", "phosphatase"
]

df = pd.read_csv(HOUSEKEEPING_BLAST_FILE, sep="\t", comment="#", header=None)
print("[INFO] Column names in `hgt_annotations.emapper.seed_orthologs`:", df.columns)

housekeeping_hits = set(df.iloc[:, 0].astype(str))
eggnog_df = pd.read_csv(ANNOTATION_FILE, sep="\t", comment="#")

if "KEGG_ko" in eggnog_df.columns:
    housekeeping_eggnog = eggnog_df[
        eggnog_df["KEGG_ko"].apply(lambda x: any(keyword in str(x).lower() for keyword in HOUSEKEEPING_KEYWORDS))
    ]["query_name"].tolist()
else:
    print("[WARNING] `KEGG_ko` column not found in annotations file. Skipping KEGG filtering.")
    housekeeping_eggnog = []

housekeeping_genes = housekeeping_hits.union(set(housekeeping_eggnog))

filtered_sequences = []
for record in SeqIO.parse(INPUT_FASTA, "fasta"):
    if record.id not in housekeeping_genes:
        filtered_sequences.append(record)

SeqIO.write(filtered_sequences, OUTPUT_FASTA, "fasta")
