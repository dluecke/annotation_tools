#!/bin/bash

# blast6togff3.sh takes blast outfmt 6 output file produces gff3 with source=blast type=match
# REQUIRES genometools gt gff3 for format check/final processing
# output gff file with same name as blast output file, but .gff extension

if [[ $# != 1 ]] || [[ ! -f $1 ]]; then
	echo "USAGE: $0 BLAST_OUTFMT6.tsv"
	echo "Outputs .gff and .tmp.gff files with otherwise same name as BLAST file to first ."
	exit
fi

OUTFILEnoext=$(basename $1 | sed 's/.tsv//')

if [ -f ${OUTFILEnoext}.gff ]; then
	echo "file ${OUTFILEnoext}.gff already exists, delete or rename and rerun"
	exit
fi

oldIFS=$IFS
IFS=$'\n'

echo '##gff-version 3' > ${OUTFILEnoext}.tmp.gff

while read LINE; do

	GENEID=$(echo $LINE | cut -f1)
	SCAFID=$(echo $LINE | cut -f2 | cut -d';' -f1) # for prolongata dovetail assembly need to remove HRSCAF from scaffold name so cut -f1 at ;
	NTERM=$(echo $LINE | cut -f9)
	CTERM=$(echo $LINE | cut -f10)
		
	if [[ $NTERM -lt $CTERM ]]; then
		START=$NTERM
		END=$CTERM
		DIR='+'
	else
		START=$CTERM
		END=$NTERM
		DIR='-'
	fi
	
	echo -e "$SCAFID\tblast\tmatch\t$START\t$END\t.\t$DIR\t.\tlocus_tag=$GENEID" >> ${OUTFILEnoext}.tmp.gff

done < $1

IFS=$oldIFS

echo
echo "Pre-formatted gff file ${OUTFILEnoext}.tmp.gff written. Should work for some browsers (eg Artemis)"
echo "Using genometools to sort and tidy into ${OUTFILEnoext}.gff"
echo 

gt gff3 -sort yes -tidy yes -o ${OUTFILEnoext}.gff ${OUTFILEnoext}.tmp.gff

