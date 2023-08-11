#!/bin/bash
# run_nucmer_serial_querys.sh runs NUCmer alignment with multiple queries to same reference
# takes -p "PREFIX", -c numeric, -r REFERENCE.fa, and -q DIRECTORY/ containing SEQ_N.fa query sequences
# runs NUCmer alignment on each SEQ_N.fa with --maxmatch, -p flag "PREFIX-SEQ_N", and given -c value
# REQUIRES MUMmer (written with version 4.0.0rc)

# read command line values
while getopts "p:c:r:q:" opt; do
	case $opt in
		p)
			PREFIX=$OPTARG
			;;
		c)
			C_VAL=$OPTARG
			;;
		r)
			REFERENCE=$OPTARG
			;;
		q)
			QUERY_DIR=$OPTARG
			;;
	esac
done

if [[ ! -n $PREFIX ]] || [[ ! -n $C_VAL ]] || [[ ! -f $REFERENCE ]] || \
   [[ ! -d $QUERY_DIR ]] || [[ ! -x $(command -v nucmer) ]]; then
   	echo "USAGE: $0 -p STRING -c NUMERIC -r REFERENCE.fa -q DIRECTORY/"
	echo "Runs series of NUCmer alignments with all DIRECTORY/*.fa files as queries onto REFERENCE.fa"
	echo "Uses NUCmer options match length -c NUMERIC and output prefix -p STRING-queryID"
	echo "REQUIRES: mummer::nucmer (written with mummer/4.0.0rc)"
	exit
fi


