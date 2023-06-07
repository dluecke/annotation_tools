#!/bin/bash

# split_assembly_by_BUSCO_dups.sh removes suspected duplicate scaffolds from assembly 
# takes RemoveList output from BUSCO_stats_by_contig.py
# writes 3 fastas:
#	ASSEMBLY-DeDup.fa, original assembly with duplicate scaffolds removed
#	ASSEMBLY-DupRemoved.fa, all removed scaffolds
#	ASSEMBLY-DupPairs.fa, retained scaffolds containing pairs for removed duplicates
# Alignment visualization between DupRemoved and DupPairs can confirm nature of duplicates
# Inspired by https://bioinformatics.stackexchange.com/a/14421 (answer on https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta)

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ $1 != *"RemoveList.tsv" ]] || [[ $2 != *".fa"* ]]; then
	echo "USAGE: $0 REMOVELIST.tsv ASSEMBLY.fa"
	echo "Reads RemoveList.tsv output from BUSCO_stats_by_contig.py"
	echo "Removes identified duplicate scaffolds from ASSEMBLY.fa"
	echo "Writes ASSEMBLY-DeDup.fa, ASSEMBLY-DupRemoved.fa, and ASSEMBLY-DupPairs.fa (for alignment visualization)"
	echo "REQUIRES: samtools"
	exit
fi


# Get strings for output filenames
ASSEMBLY="${2%.*}" # removes file extension
OUTFILE_DEDUP=$ASSEMBLY"-DeDup.fa"
OUTFILE_DUPREMOVE=$ASSEMBLY"-DupRemoved.fa"
OUTFILE_DUPPAIRS=$ASSEMBLY"-DupPairs.fa"


# Ensure have a samtools fasta index
if [[ ! -f $2".fai" ]]; then
	samtools faidx $2
fi


# Make arrays for Scaffold IDs from RemoveList.tsv columns 

# remove column already in order and one per line
remove_ids=($(cut -f1 $1 | tail -n +2 ))

# keep (pairs) column has comma separated lists, and want to order by scaffold number 
pairs_ids=($(cut -f2 $1 | sed 's/, /\n/' | \
	awk 'NR > 1 {print length, $0 }' | \
	sort -n -k1,1 -k2 | cut -d' ' -f2 | uniq ))

# for list of all scaffolds to keep need to use grep -v 
keep_ids=($(awk '{print $1}' $2.fai | \
	grep -v -F "$(printf '%s\n' "${remove_ids[@]}")" ))


# Use samtools faidx to extract relevant scaffolds for each and write to new files
samtools faidx -o $OUTFILE_DEDUP $2 "${keep_ids[@]}"
samtools faidx -o $OUTFILE_DUPREMOVE $2 "${remove_ids[@]}"
samtools faidx -o $OUTFILE_DUPPAIRS $2 "${pairs_ids[@]}"

