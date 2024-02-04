#!/bin/bash

# annotation_BLAST_scoring.sh shell to run blast6togff3.sh, merge_gff3.sh, annotation_vs_BLASTregions.sh
# $1 assembly.fasta - assembly for annotations
# $2 annotations.fofn - file of annotation file relative paths
# $3 proteins.fasta - protein sequences to score annotation's BLAST hit coverage
# $4 ProtIDstring - shorter string to identify protein set

# REQUIRES genometools gt gff3, bedtools merge

if [[ $# != 4 ]] || [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ $2 != *".fofn" ]] || [[ ! -f $3 ]]; then
	echo "USAGE: $0 ASSEMBLY.fa ANNOTATIONS.fofn PROTEINS.fa Prot_ID_string"
	echo "checks if regions with BLAST hits to genes of interest are covered by annotation"
	echo "REQUIRES: genometools gt gff3, bedtools merge"
	exit
fi

ASSEMBLYNAME=$(basename $1 | cut -d'.' -f1)
echo "ASSEMBLYNAME:"
echo $ASSEMBLYNAME

### BLAST PROTEINS TO ASSEMBLY ###

# make blast database
makeblastdb -dbtype nucl -in $1

BLASTOUTNAME=$4-to-$ASSEMBLYNAME-e001.tsv
tblastn -db $1 -query $3 -evalue 0.001 -outfmt 6 -out $BLASTOUTNAME

echo "BLAST results finished: $BLASTOUTNAME"
wc $BLASTOUTNAME
echo

### CONVERT BLAST OUTPUT TO GFF ###

# make long gff - multiple entries per region

~/scripts/blast6togff3.sh $BLASTOUTNAME

# long gff filename based on BLASTOUTNAME
BLASTnoext=$(echo $BLASTOUTNAME | sed 's/.tsv//')
LONG_GFF=$BLASTnoext.gff

# make short gff - regions unique

~/scripts/merge_gff3.sh $LONG_GFF

# long gff filename based on BLASTOUTNAME
GFF=$BLASTnoext.merge.gff


### COMPARE ALL ANNOTATIONS TO BLAST GFF ###
oldIFS=$IFS
IFS=$'\n'
while read annotation; do

	~/scripts/annotation_vs_BLASTregions.sh $GFF $annotation $4
	
done < $2
IFS=$oldIFS



