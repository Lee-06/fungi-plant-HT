#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Calculate HGT Index by comparing Fungi vs Plant Bitscores.")
parser.add_argument("--fungi_results", required=True, help="Filtered BLAST results against Fungi (from Script 2)")
parser.add_argument("--plant_results", required=True, help="BLAST results against Plants (from Script 4)")
parser.add_argument("--output", default="fungi_vs_plant_comparison.tsv", help="Output comparison file")
parser.add_argument("--candidates_out", default="hgt_candidates.tsv", help="Final HGT candidates file")

args = parser.parse_args()

print("[INFO] Loading Fungi results...")
fungi_df = pd.read_csv(args.fungi_results, sep="\t")

print("[INFO] Loading Plant results...")
plant_df = pd.read_csv(args.plant_results, sep="\t")

fungi_best = fungi_df.sort_values("bitscore", ascending=False).groupby("qseqid", as_index=False).first()
fungi_best = fungi_best.rename(columns=lambda x: x + "_fungi" if x != "qseqid" else x)

plant_best = plant_df.sort_values("bitscore", ascending=False).groupby("qseqid", as_index=False).first()
plant_best = plant_best.rename(columns=lambda x: x + "_plant" if x != "qseqid" else x)

print("[INFO] Merging and calculating HT Index...")
comparison = pd.merge(fungi_best, plant_best, on="qseqid", how="outer")

def calculate_metrics(row):
    bf = row.get("bitscore_fungi", 0)
    bp = row.get("bitscore_plant", 0)
    if pd.isna(bf): bf = 0
    if pd.isna(bp): bp = 0
    h_index = bf - bp
    ratio = bf / bp if bp > 0 else bf 
    return pd.Series([h_index, ratio], index=['h_index', 'score_ratio'])

comparison[["h_index", "score_ratio"]] = comparison.apply(calculate_metrics, axis=1)
comparison_sorted = comparison.sort_values("h_index", ascending=False)

comparison_sorted.to_csv(args.output, sep="\t", index=False)

candidates = comparison_sorted[comparison_sorted["h_index"] > 0]
candidates.to_csv(args.candidates_out, sep="\t", index=False)

print(f"[SUCCESS] Comparison saved to {args.output}")
print(f"[INFO] {len(candidates)} potential candidates (h_index > 0) saved to {args.candidates_out}")
