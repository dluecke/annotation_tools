#!/bin/bash

# annotation_vs_BLASTregions.sh scores annotation based on blast hit regions
# REQUIRES bedtools intersect (Version: v2.29.2 when written)
# $1: reference gff of UNIQUE regions with blast hits (blast6togff3.sh and merge_gff3.sh) 
# $2: annotation file to be scored
# (opt) $3: description of blast hits for output filename (otherwise use blast gff filename to first .)

if [[ ! -f $1 ]] || [[ ! -f $2 ]]; then
	echo "USAGE: $0 MERGED_BLAST_REGIONS.gff ANNOTATION_TO_SCORE.gff BLAST_ID_STRING"
	echo "Outputs <ANNOTATION_TO_SCORE.gff>_vs_<BLAST_ID_STRING>.tsv with line for each entry in BLAST gff and largest overlap found"
	exit
fi

if [[ ! -z $3 ]]; then
	BLAST_ID=$3
else
	BLAST_ID=$(basename $1 | cut -d'.' -f1)
fi

ANNOTATIONnoext=$(basename $2)

OUTFILE="${ANNOTATIONnoext}-vs-${BLAST_ID}.tsv"

if [ -f $OUTFILE ]; then
	echo "file $OUTFILE already exists, delete or rename and rerun"
	exit
fi

# Call bedtools intersect on gff without contig sections (these give full hits but aren't what we want)
# Take relevant columns of output (position of BLAST hit, annotation type, overlap position and length)
# One line per BLAST region (annotation feature with longest overlap), direct to OUTFILE
bedtools intersect -a $1 -b <(sed '/[[:space:]]contig[[:space:]]/d' $2) -wao \
	| cut -f1,3-5,7,11-14,19 | sort -rn -k10,10 | awk -F"\t" '!_[$3]++' | sort -k1,1 -k3,3n > $OUTFILE

FOUND_COUNT=$(cut -f10 $OUTFILE | grep "^[^0]" -c)
TOTAL_COUNT=$(wc -l < $OUTFILE | tr -d ' ')

echo
echo "$FOUND_COUNT out of $TOTAL_COUNT regions with BLAST hits are annotated"
echo "See $OUTFILE for details of each region"
echo

