#!/bin/bash

# convert_gnuplot_to_tsv.sh takes gnuplot data file and converts to tsv
# intended to take mummerplot data files and import to R for custom plotting
# mummerplot writes data files with extensions .fplot and .rplot (forward and reverse)
# gnuplot data file format:
#  (blank line)
#  R1 Q1 similarity
#  R2 Q2 similarity
#  (blank line)
# writes tsv NAME.[fr]plot.tsv with format:
#  	R1	Q1	R2	Q2
# drops similarity score

if [[ ! -f $1 ]]; then
	echo "USAGE: $0 plotdata.[fr]plot"
	echo "Takes gnuplot data files from mummerplot either .fplot or .rplot"
	echo "Converts to tsv file for R plotting"
	echo "Writes file plotdata.[fr]plot.tsv"
	exit
fi

TSVFILE=$1".tsv"

if [[ -f $TSVFILE ]]; then
	echo "$TSVFILE exists, delete or rename then rerun"
	exit
fi

echo -e "R1\tQ1\tR2\tQ2" > $TSVFILE

awk 'NR > 3 {
	if (NF > 0) {
		if (line != "") {
			line = line "\t" $1 "\t" $2; 
		} else {
			line = $1 "\t" $2; 
		}
	} 
	if(NF == 0) { 
		if (line != "") { 
			print line; line = ""; 
		} 
	}
}
END {
	if (line != "") {
		print line; 
	}
}'
		
		
		
		
		
		