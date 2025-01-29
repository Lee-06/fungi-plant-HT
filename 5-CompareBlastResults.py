#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lee Mariault
"""
import pandas as pd
import numpy as np

FUNGI_RESULTS_FILE     = "filtered_blast_results_with_fungi.tsv"  
PLANT_ALIGNMENT_FILE   = "plant_alignment_results.tsv"            
FUNGI_VS_PLANT_FILE    = "fungi_vs_plant_comparison.tsv"          

columns_fungi = [
    "qseqid", "sseqid", "pident", "length", 
    "mismatch", "gapopen", "qstart", "qend", 
    "sstart", "send", "evalue", "bitscore", 
    "fungi_genome"
]

fungi_df = pd.read_csv(FUNGI_RESULTS_FILE, sep="\t", names=columns_fungi, header=0, engine="python")

columns_plant = [
    "qseqid", "sseqid", "pident", "length", 
    "evalue", "bitscore", "plant_genome"
]
plant_df = pd.read_csv(PLANT_ALIGNMENT_FILE, sep="\t", names=columns_plant, header=0, engine="python")

fungi_df_sorted = fungi_df.sort_values("bitscore", ascending=False)
best_fungal_hits = fungi_df_sorted.groupby("qseqid", as_index=False).first()

best_fungal_hits = best_fungal_hits.rename(columns={
    "sseqid":   "sseqid_fungi",
    "pident":   "pident_fungi",
    "length":   "length_fungi",
    "mismatch": "mismatch_fungi",
    "gapopen":  "gapopen_fungi",
    "qstart":   "qstart_fungi",
    "qend":     "qend_fungi",
    "sstart":   "sstart_fungi",
    "send":     "send_fungi",
    "evalue":   "evalue_fungi",
    "bitscore": "bitscore_fungi",
    "fungi_genome": "fungi_genome"
})

plant_df_sorted = plant_df.sort_values("bitscore", ascending=False)
best_plant_hits = plant_df_sorted.groupby("qseqid", as_index=False).first()

best_plant_hits = best_plant_hits.rename(columns={
    "sseqid":   "sseqid_plant",
    "pident":   "pident_plant",
    "length":   "length_plant",
    "evalue":   "evalue_plant",
    "bitscore": "bitscore_plant",
    "plant_genome": "plant_genome"
})

comparison = pd.merge(
    best_fungal_hits,
    best_plant_hits,
    on="qseqid",
    how="outer"  
)

def kingdom_predominance(row):
    bf = row.get("bitscore_fungi", np.nan)
    bp = row.get("bitscore_plant", np.nan)
    if pd.isna(bf) and pd.isna(bp):
        return "tie"  
    elif pd.isna(bf) and not pd.isna(bp):
        return "plant"  
    elif pd.isna(bp) and not pd.isna(bf):
        return "fungi"  
    else:
        if bf > bp:
            return "fungi"
        elif bp > bf:
            return "plant"
        else:
            return "tie"

comparison["predominance"] = comparison.apply(kingdom_predominance, axis=1)
comparison.to_csv(FUNGI_VS_PLANT_FILE, sep="\t", index=False)
print(f"[INFO] Final comparison saved to '{FUNGI_VS_PLANT_FILE}'.")
