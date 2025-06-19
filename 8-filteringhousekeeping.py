#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault

"""

import sys
import pandas as pd
from Bio import SeqIO

# ── User paths ────────────────────────────────────────────────────────────────
ANNOTATION_FILE         = "hgt_annotations.emapper.annotations"
BLAST_FILE              = "hgt_annotations.emapper.seed_orthologs"
INPUT_FASTA             = "hgt_clusters.fasta"
OUTPUT_FASTA            = "hgt_filtered.fasta"

# ── Housekeeping keywords ─────────────────────────────────────────────────────
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
# lowercase for case-insensitive matching
keywords_lower = [kw.lower() for kw in HOUSEKEEPING_KEYWORDS]

try:
    blast_df = pd.read_csv(BLAST_FILE, sep="\t", header=None,
                           comment="#", dtype=str)
except FileNotFoundError:
    sys.exit(f"[ERROR] Cannot find seed-orthologs file: {BLAST_FILE!r}")
print("[INFO] Columns in seed-ortholog file:", blast_df.columns.tolist())
all_hits = set(blast_df.iloc[:, 0].astype(str))

header_idx = None
header_line = None
with open(ANNOTATION_FILE, 'r') as af:
    for idx, line in enumerate(af):
        if line.startswith("#query_name"):
            header_line = line.lstrip("#").strip()
            header_idx = idx
            break

if header_idx is None:
    sys.exit(f"[ERROR] Could not find a line beginning `#query_name` in {ANNOTATION_FILE}")

column_names = header_line.split("\t")

try:
    eggnog_df = pd.read_csv(
        ANNOTATION_FILE,
        sep="\t",
        skiprows=header_idx + 1,
        header=None,
        names=column_names,
        dtype=str
    )
except Exception as e:
    sys.exit(f"[ERROR] Failed to parse annotations: {e}")

print("[INFO] Columns in annotation file:", eggnog_df.columns.tolist())

desc_col = next((c for c in eggnog_df.columns if "free text" in c.lower()), None)
to_scan = ["Preferred_name"]
if desc_col:
    to_scan.append(desc_col)

hk_candidates = []
for _, row in eggnog_df.iterrows():
    combined = []
    for col in to_scan:
        val = row.get(col, "")
        if pd.notna(val):
            combined.append(str(val))
    text = " ".join(combined).lower()
    if any(kw in text for kw in keywords_lower):
        hk_candidates.append(str(row["query_name"]))

print(f"[INFO] {len(hk_candidates)} queries match housekeeping keywords in name/desc")

housekeeping_genes = set(hk_candidates).intersection(all_hits)
print(f"[INFO] Will remove {len(housekeeping_genes)} housekeeping genes total")

kept = []
for rec in SeqIO.parse(INPUT_FASTA, "fasta"):
    if rec.id not in housekeeping_genes:
        kept.append(rec)

SeqIO.write(kept, OUTPUT_FASTA, "fasta")
print(f"[INFO] Removed {len(housekeeping_genes)} sequences; saved {len(kept)} in {OUTPUT_FASTA}")
