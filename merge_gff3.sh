#!/bin/bash

# merge_gff3.sh takes gff with many overlapping regions (eg BLAST hits)
# merges overlapping entries, including within 500bp (-d option in bedtools merge), into single entry
# outputs gff3 file formatted as input file (sequence-regions defined at start)
# REQUIRES bedtools merge 
# expects gff3 output from genometools gff3 -tidy -sort, as produced by blast6t0gff3.sh

if [[ $# != 1 ]] || [[ ! -f $1 ]]; then
	echo "USAGE: $0 BLASTtoGFF3_OUTPUT.gff3"
	echo "Outputs .merge.gff and .tmp.bam files with otherwise same name as BLAST file to first ."
	exit
fi

OUTFILEnoext=$(basename $1 | sed 's/.gff//')

if [ -f ${OUTFILEnoext}.merge.gff ]; then
	echo "file ${OUTFILEnoext}.merge.gff already exists, delete or rename and rerun"
	exit
fi

echo
echo "using bedtools merge to produce ${OUTFILEnoext}.tmp.bam file for merged regions with scaffold, start, end, strand, and locusIDs"

bedtools merge -d 500 -c 7,9 -o distinct,distinct -s -i $1 > ${OUTFILEnoext}.tmp.bam

echo
echo "start writing ${OUTFILEnoext}.merge.gff with initial commented block from $1 (gff version and sequence-regions)"

FirstEntryLine=$(grep -n "^[^#]" $1 | head -n1 | cut -d':' -f1)

head -n$((FirstEntryLine-1)) $1 > ${OUTFILEnoext}.merge.gff

echo
echo "formatting ${OUTFILEnoext}.tmp.bam into gff3 format with source=blast type=match, appending to ${OUTFILEnoext}.merge.gff"

# add 1 to positions for 0 to 1 indexing
awk 'FS=OFS="\t" { print $1, "match", "blast", $2+1, $3+1, ".", $4, ".", $5 }' ${OUTFILEnoext}.tmp.bam >> ${OUTFILEnoext}.merge.gff

