#!/bin/bash
# get_gene_regions.sh takes ASSEMBLY.fa and ANNOTATION.gff
# extracts all oriented genome regions with "gene" annotation feature
# writes new fasta ASSEMBLY-GeneRegions.fa in same directory as ASSEMBLY, 
#  which contains all regions with annotated genes in correct orientation
#  sequences labeled with annotation ID (eg "ID=geneN"),
#  corresponding scaffold and coordinates, and orientation compared to assembly
# also writes new gff ANNOTATION-GeneRegions.gff in same dir as ANNOTATION,
#  which contains only lines with "gene" feature, but replaces 3rd column ("gene") with
#  first field of 9th column (eg ID=geneN) 
#   - enables labeling of gene ID in each region fasta header
# REQUIRES bedtools2 (written with version 2.29.2)

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ ! -x $(command -v bedtools getfasta) ]]; then
	echo "USAGE: $0 ASSEMBLY.fa ANNOTATION.gff"
	echo "Writes ASSEMBLY-GeneRegions.fa with all annotated gene regions in correct orientation"
	echo "Outputs to same directory as original assembly"
	echo "REQUIRES: bedtools getfasta"
	exit
fi

# get relevant file names and directory paths from input files
FASTAFILE=$1
FASTADIR=$(dirname $(readlink -f $FASTAFILE))
FASTAFILENAME=$(basename $FASTAFILE)
FASTAFILENOEXT=$(echo "$FASTAFILENAME" | rev | cut --complement -f1 -d'.' | rev)
FASTAFILEABS=$FASTADIR/$FASTAFILENAME
FASTAOUT=$FASTAFILENOEXT"-GeneRegions.fa"

GFFFILE=$2
GFFDIR=$(dirname $(readlink -f $GFFFILE))
GFFFILENAME=$(basename $GFFFILE)
GFFFILENOEXT=$(echo "$GFFFILENAME" | rev | cut --complement -f1 -d'.' | rev)
GFFFILEABS=$GFFDIR/$GFFFILENAME
GFFOUT=$GFFFILENOEXT"-GeneRegions.gff"

# write new GFF with only gene feature lines, replacing 3rd column with geneID (first item in column 9)
cd $GFFDIR
awk '$3== "gene" { split($9, id, ";"); 
	print $1 "\t" $2 "\t" id[1] "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' \
	$GFFFILE > $GFFOUT

# write gene region fasta using -s option in getfasta for strandedness 
cd $FASTADIR
bedtools getfasta -fi $FASTAFILE -bed $GFFFILEABS -s -name -fo $FASTAOUT



