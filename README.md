__**Overview**__

This repository contains the bioinformatics pipeline used to identify, filter, and validate Horizontal Transfer (HT) candidates between fungi and plants. The pipeline performs a massive all-vs-all genomic comparison, filters for high-confidence transfer events, removes contaminants/housekeeping genes, and validates candidates using phylogenetic reconstruction.

In order for the pipeline to function correctly, it is highly recommended to have a dataset containing genomes from a wide taxonomical range (at least 1 species per family, dozens of families represented per orders).

__**System Requirements**__

*Dependencies*
The pipeline requires the following tools to be installed and in your system $PATH:
BLAST+ (specifically blastn, makeblastdb)
hs-blastn (Highly recommended for the initial genome-wide search)
CD-HIT (for clustering)
EggNOG-mapper (for functional annotation)
MAFFT (for multiple sequence alignment)
IQ-TREE (for phylogenetic tree inference)
Diamond (used by EggNOG-mapper)

*Python Libraries*
pandas
biopython
numpy

__**Usage Guide**__

*Step 1: Genome-Wide Homology Search*
./1-BlastWholeGenomes.sh ./data/plant_genomes ./data/fungi_genomes

*Step 2: Primary Filtering*
python 2-filter_blast_results.py \
    --blast_dir ./blastresults \
    --fungi_fai ./data/fungi_indices \
    --plant_fai ./data/plant_indices
    
*Step 3: Extract & Check Distribution*
python 3-extractfasta.py

python 4-FindNonUbiquitousSequences.py \
    -s ./selected_sequences \
    -p ./data/plant_genomes \
    -j 10 -t 4

*Step 4: Calculate HT Index*
python 5-CompareBlastResults.py

*Step 5: Functional Annotation & Cleaning*
python 6-extractHTcandidates.py

./7-cluster_and_annotate_candidates.sh /path/to/eggnog_database

python 8-filteringhousekeeping.py

*Step 6: Phylogenetic Validation*
python 9-build_phylogenies.py \
    --input hgt_filtered.fasta \
    --database /path/to/ncbi_nt_db \
    --outdir ./phylogenies

__**Output**__

hgt_candidates.tsv: Table of potential HT events with scores.
phylogenies/: Directory containing .treefile (Newick trees) and alignments for every candidate. In the current version these trees must be inspected for fungal nesting within plant clades (or vice-versa). *(In future versions this verification step is automated by RANGER-DTL, a tree reconciliation software)*
