#!/bin/bash

#minimum identity
minid=80
#minimum alignment length
minln=500
#minimum scaffold length
minsf=2000

#Separating blast results into distinct files for each fungi species
for blastfile in *.blast;do
	cut -f1 -d'#' $blastfile | sort -u > fungilist.txt
	while read -r line;do
		grep $line $blastfile | sed 's/^.*#//' > 1-separated/${line}\#${blastfile#*.fasta_};
	done < fungilist.txt
done

#Creating summary file from blast results that also include scaffold length
for startingletter in {A..Z};do
	for separatedblast in 1-separated/${startingletter}*.blast; do
		funginame=$(echo "$(basename "$separatedblast")" | cut -d"#" -f1)
		plantname=$(echo "$(basename "$separatedblast" .blast)" | cut -d"#" -f2)
		fungilength=$(echo "/home/lee/Bureau/fungi_genomes/${funginame}*_AssemblyScaffolds.fasta.length")
		plantlength=$(echo "/home/lee/Bureau/Genome_Plant_All/${plantname}.fasta.length")
		while IFS=$'\t' read -r fungiscaffold plantscaffold identity alnlength mismatch gaps qstart qend sstart send e_value bitscore
		do
			fungiscaffoldname=${fungiscaffold%:*}
			fungiscaffoldinterval=${fungiscaffold#*:}
			fungiscaffoldstartpos=${fungiscaffoldinterval%-*}
			fungscafflength=$(grep -w -m 1 $fungiscaffoldname $fungilength | cut -d' ' -f1)
			plantscafflength=$(grep -w -m 1 $plantscaffold $plantlength | cut -d' ' -f1)
			fungalnrealstart=$((fungiscaffoldstartpos+qstart))
			fungalnrealend=$((fungiscaffoldstartpos+qend))
			str=$funginame" "$plantname" "$identity" "$alnlength" "$fungiscaffoldname" "$fungalnrealstart" "$fungalnrealend" "$fungscafflength" "$plantscaffold" "$sstart" "$send" "$plantscafflength" "$e_value" "$bitscore
			echo $str >> 2-summary/${startingletter}_fungi_VS_All_Plants.blast
		done < "$separatedblast"
	done
done

#Fixing some blast results that have plant scaffolds alignment coordinates backwards
for summedblast in 2-summary/*.blast;do
	awk '$10 > $11 { temp = $11; $11 = $10; $10 = temp } 1' $summedblast | nl -s" " -n rz -w 9 > 3-reordered/reordered_${summedblast}
done

#Applying my filters to keep high identity, minimum length alignments on scaffolds with a minimum length as well
for reorderedblast in 3-reordered/*.blast; do
	reorderedfilename=$(echo "$(basename "$reorderedblast")")
	awk -v i="$minid" -v l="$minln" -v s="$minsf" '$4>i && $5>=l && $9>s && $13>s' $file | sed 's/ /\t/g' > 4-filtered/${minid}id.aln${minln}bp.scafsup${minsf}bp.${reorderedfilename}.fungi.blast
done

#Merging and simplifying blast results into BED format
for filteredblast in 4-filtered/*.blast; do
	filteredfilename=$(echo "$(basename "$filteredblast")")
	sed 's/ /\t/g' $filteredblast | sort -u | awk 'OFS="\t"{print $6,$7,$8,$2}' | sort -k1,1 -k2,2n | mergeBed -d 1 -c 4 -o distinct -i - > 5-merged/${filteredfilename}.merged.bed
done



########## BELOW IS A WORK IN PROGRESS ###############



cut -f4 5-merged/*.merged.bed | sort -u > 5-merged/FungiList.txt

while IFS=$'\t' read -r species; do
	grep $species 5-merged/*.merged.bed > 6-smallbeds/${species%fasta}small.bed
done < 5-merged/FungiList.txt

for fastafile in mergedbedfiles/EpichloeVSPlants/*.small.bed; do
	actualfilename=$(echo "$(basename "$fastafile")")
	fastaFromBed -fi /run/user/1000/gvfs/ftp:host=nas-lgdp.local/Eq_Panaud/Moaine/fungi_genomes/"${actualfilename%.small.bed}.fasta" -bed $fastafile -fo ${fastafile%.small.bed}.PlantAligned.fasta
done

mkdir mergedbedfiles/EpichloeVSPlants/merged

mv mergedbedfiles/EpichloeVSPlants/*.PlantAligned.fasta .

for scaffoldfile in *.PlantAligned.fasta; do
	awk '/>/{sub(">","&"FILENAME"#");sub(/\.PlantAligned.fasta/,x)}1' $scaffoldfile >> mergedbedfiles/EpichloeVSPlants/merged/AllEpichloeAlignedScaffolds.fasta
done

mv *.PlantAligned.fasta mergedbedfiles/EpichloeVSPlants/

mkdir mergedbedfiles/EpichloeVSPlants/merged/silvablastresults

blastn -db /home/lee/Bureau/sequencefiltering/DB/silvaribo.fasta -query mergedbedfiles/EpichloeVSPlants/merged/AllEpichloeAlignedScaffolds.fasta -outfmt 6 -out mergedbedfiles/EpichloeVSPlants/merged/silvablastresults/AllEpichloeSILVA.blast -evalue 1e-20 -num_threads 20

cut -f1 mergedbedfiles/EpichloeVSPlants/merged/silvablastresults/AllEpichloeSILVA.blast > mergedbedfiles/EpichloeVSPlants/merged/silvablastresults/silvariboscaffolds.txt

python /home/lee/Bureau/sequencefiltering/filter_fasta_by_list_of_headers.py mergedbedfiles/EpichloeVSPlants/merged/AllEpichloeAlignedScaffolds.fasta mergedbedfiles/EpichloeVSPlants/merged/silvablastresults/silvariboscaffolds.txt > mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta

makeblastdb -in mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta -dbtype nucl
blastn -db mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta -query mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta -out mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta.blast -outfmt 6 -evalue 1e-20 -num_threads 20
silix -i 0.4 -r 0.4 mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta.blast -f FAM > mergedbedfiles/EpichloeVSPlants/merged/EpichloeClusters.txt
fastareformat mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.fasta > mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.reformated.fasta

mkdir mergedbedfiles/EpichloeVSPlants/merged/clustering

silix-split mergedbedfiles/EpichloeVSPlants/merged/silvafiltered.AllEpichloeAlignedScaffolds.reformated.fasta mergedbedfiles/EpichloeVSPlants/merged/EpichloeClusters.txt -o mergedbedfiles/EpichloeVSPlants/merged/clustering/ -n 1

for i in mergedbedfiles/EpichloeVSPlants/merged/clustering/*.fasta; do 
	cd-hit -i $i -o $i.cdhit -c 0.7 
done

cat mergedbedfiles/EpichloeVSPlants/merged/clustering/*.fasta > silvafiltered.AllEpichloeAlignedScaffolds.reformated_ALLFAMILIES.fasta
