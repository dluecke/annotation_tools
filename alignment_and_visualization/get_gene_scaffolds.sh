#!/bin/bash
# get_gene_scaffolds.sh takes ASSEMBLY.fa and ANNOTATION.gff
# writes new fasta ASSEMBLY-GeneScafs.fa, which
# contains all scaffolds with annotated genes
# REQUIRES samtools (written with version 1.15.1)

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -x $(command -v samtools faidx) ]]; then
	echo "USAGE: $0 ASSEMBLY.fa ANNOTATION.gff"
	echo "Writes assembly subset ASSEMBLY-GeneScafs.fa (all scaffolds with annotated genes)"
	echo "Outputs to same directory as original assembly"
	echo "REQUIRES: samtools faidx"
	exit
fi

# get relevant file and directory names from input file
INFILE=$1
INDIR=$(dirname $INFILE)
FILENAME=$(basename $INFILE)
NOEXTNAME=$(echo "$FILENAME" | rev | cut --complement -f1 -d'.' | rev)

# get sequences with gene annotations from provided gff
# sorted by scaffold number assuming Scaffold_N format
cut -f1-3 $2 | grep -w gene | cut -f1 | sort -u | \
     awk '{print length, $0 }' | sort -k1,1 | \
     cut -f2 -d' ' > $INDIR/$NOEXTNAME-GeneScafNames.txt

# go to assembly directory for samtools output
cd $INDIR
samtools faidx -o $NOEXTNAME-GeneScafs.fa --region $NOEXTNAME-GeneScafNames.txt $FILENAME
