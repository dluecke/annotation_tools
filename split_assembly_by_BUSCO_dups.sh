#!/bin/bash

# split_assembly_by_BUSCO_dups.sh removes suspected duplicate scaffolds from assembly 
# takes RemoveList output from BUSCO_stats_by_contig.py
# writes 3 fastas:
#	ASSEMBLY-DeDup.fa, original assembly with duplicate scaffolds removed
#	ASSEMBLY-DupRemoved.fa, all removed scaffolds
#	ASSEMBLY-DupPairs.fa, retained scaffolds containing pairs for removed duplicates
# Alignment visualization between DupRemoved and DupPairs can confirm nature of duplicates
# Inspired by https://bioinformatics.stackexchange.com/a/14421 (answer on https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta)

if [[ ! -f $1 ]] || [[ ! -f $2 ]] || [[ $1 != *"RemoveList"*".tsv" ]] || [[ $2 != *".fa"* ]]; then
	echo "USAGE: $0 REMOVELIST.tsv ASSEMBLY.fa [ANNOTATION.gff]"
	echo "Reads RemoveList.tsv output from BUSCO_stats_by_contig.py"
	echo "Removes identified duplicate scaffolds from ASSEMBLY.fa"
	echo "Writes ASSEMBLY-DeDup.fa, ASSEMBLY-DupRemoved.fa, and ASSEMBLY-DupPairs.fa (for alignment visualization)"
	echo "If ANNOTATION.gff included also writes ANNOTATION-DeDup.gff and ANNOTATION-DupRemoved.gff"
	echo "REQUIRES: samtools"
	exit
fi


# Get strings for output filenames
ASSEMBLY="${2%.*}" # removes file extension
OUTFILE_DEDUP=$ASSEMBLY"-DeDup.fa"
OUTFILE_DUPREMOVE=$ASSEMBLY"-DupRemoved.fa"
OUTFILE_DUPPAIRS=$ASSEMBLY"-DupPairs.fa"

# Make annotation output filename if file included (DupPairs included for alignments don't need annotation)
if [[ -f $3 ]]; then
	ANNOTATION="${3%.*}"
	ANNOTATION_DEDUP=$ANNOTATION"-DeDup.gff"
	ANNOTATION_DUPREMOVE=$ANNOTATION"-DupRemoved.gff"
fi

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
	
# make a sed command to delete exact word matches for remove ids
SED_DELETE_CMD=""
for pattern in "${remove_ids[@]}"; do
	SED_DELETE_CMD+=" /${pattern}\\b/d;"
done

# make a sed command to print exact word matches for remove ids
SED_PRINT_CMD=""
for pattern in "${remove_ids[@]}"; do
	SED_PRINT_CMD+=" /${pattern}\\b/p;"
done

# use sed delete command to remove appropriate scaffolds from full scaffold list pulled from fasta index
keep_ids=($(awk '{print $1}' $2.fai | sed "$SED_DELETE_CMD" ))


# Use samtools faidx to extract relevant scaffolds for each and write to new files
samtools faidx -o $OUTFILE_DEDUP $2 "${keep_ids[@]}"
samtools faidx -o $OUTFILE_DUPREMOVE $2 "${remove_ids[@]}"
samtools faidx -o $OUTFILE_DUPPAIRS $2 "${pairs_ids[@]}"

# also split annotation if included
if [[ -f $3 ]]; then
	sed "$SED_DELETE_CMD" $3 > $ANNOTATION_DEDUP
	sed "$SED_PRINT_CMD" $3 > $ANNOTATION_DUPREMOVE
fi
