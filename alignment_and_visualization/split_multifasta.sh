#!/bin/bash
# split_multifasta.sh takes multifasta file and writes separate .fa for each sequence
# writes new files in directory NAME/ based on input fasta NAME.fa[sta]
# REQUIRES samtools (written with version 1.15.1)

if [[ ! -f $1 ]] || [[ ! -x $(command -v samtools faidx) ]]; then
	echo "USAGE: $0 MULTIFASTA.fa"
	echo "Writes new .fa file for each sequence in MULTIFASTA.fa"
	echo "Makes new directory MULTIFASTA/ with new .fa files"
	echo "New files named MULTIFASTA-seqname.fa"
	echo "REQUIRES: samtools faidx"
	exit
fi

# get relevant file and directory names from input file
INFILE=$1
INDIR=$(dirname $INFILE)
FILENAME=$(basename $INFILE)
NOEXTNAME=$(echo "$FILENAME" | rev | cut --complement -f1 -d'.' | rev)
NEWDIR="$INDIR/$NOEXTNAME"

# make new subdirectory and go there (will make samtools command simpler)
mkdir $NEWDIR
cd $NEWDIR

# get names of sequences from multifasta
grep ">" ../$FILENAME | tr -d ">" > SeqNames.txt

# loop through seqnames, use samtools faidx to extract each sequence into new fasta file
while read SEQ; do
	samtools faidx -o $NOEXTNAME"-"$SEQ".fa" ../$FILENAME $SEQ
done < SeqNames.txt
