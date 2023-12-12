#!/bin/bash

# convert_gnuplot_to_tsv.sh takes gnuplot data files and converts to tsv
# intended to take mummerplot data files and import to R for custom plotting
# R script alignment_region_dotplot.R takes these outputs

# mummerplot writes data files with extensions .fplot and .rplot (forward and reverse)
# gnuplot data file format:
#  (blank line)
#  R1 Q1 similarity
#  R2 Q2 similarity
#  (blank line)
# writes tsv NAME.plot.tsv with format:
#  	R1	Q1	R2	Q2
# drops similarity score

# also finds mummerplot gp plot file "NAME.gp" and extracts scaffold breakpoints
# outputs NAME.breaks.tsv with gp script line, seq name, and start position in alignment
# need for R to mark scaffold demarcations

if [[ ! -f $1.fplot ]]; then
	echo "USAGE: $0 FILENAMEnoext"
	echo "Reads gnuplot data files from mummerplot FILENAME.fplot and FILENAME.rplot"
	echo "Converts to tsv file for R plotting"
	echo "Writes file plotdata.plot.tsv"
	echo "Also reads FILENAME.gp and writes FILENAME.breaks.tsv for R plotting"
	exit
fi

TSVFILE=$1".plot.tsv"
BREAKSFILE=$1".breaks.tsv"

if [[ -f $TSVFILE ]]; then
	echo "$TSVFILE exists, delete or rename then rerun"
	exit
fi

# headers for tsv
echo -e "R1\tQ1\tR2\tQ2" > $TSVFILE

# append forward plot data
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
}' $1.fplot >> $TSVFILE

# append reverse plot data
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
}' $1.rplot >> $TSVFILE

# get the scaffold breakpoints from the .gp file if file can be found
if [[ -f $1.gp ]]; then
	echo -e "gpLine\tSeqName\tPosition" > $BREAKSFILE
	# lines with scaffold breakpoints start with space then " character
	# have to remove 2 extraneous '^ "' lines at end
	grep -n '^ "' $1.gp | sed '$d' | sed '$d' | \
		tr -d ':' | tr -d ',' | \
		awk '{print $1 "\t" $2 "\t" $3}' >> $BREAKSFILE
fi
	


