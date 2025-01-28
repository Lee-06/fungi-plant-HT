#!/bin/bash
for plantgenome in [plant_genomes]/*.fasta;do
	plantfilename=$(basename "$plantgenome")
	for fungigenome in [fungi_genomes]/*.fasta;do
		fungifilename=$(basename "$fungigenome")
		blastn -db $plantgenome -query $fungigenome -num_threads 20 -evalue 1e-20 -outfmt 6 -out blastresults/${plantfilename}_VS_${fungifilename}.blast
	done
 done
