#!/usr/bin/env python3

# remove_existing_gff_features.py writes gff without genes found by gffcompare reference comparison
# takes ($1) NEW.gff annotation and ($2) gffcmp.NEW.gff.refmap 
# probably want to strip NEW.gff of non-feature lines (written using Liftoff output)
# writes lines from NEW.gff without match in refmap qry_id_list (specifically field before | )
# outputs NEW-unique.gff (these lines can be added to reference gff without introducing redundancy)
# use gffcompare -r REFERENCE.gff NEW.gff to get refmap file

import sys
import pandas as pd
import csv

GFF = pd.read_csv(sys.argv[1], sep='\t', header=None)

REFMAP = pd.read_csv(sys.argv[2], sep='\t')

OUTGFF = sys.argv[1] + "-unique.gff"


# Need list of gene IDs from files
# GFF gene=...; substring extracted from 9th column
GFFgeneIDs = GFF[8].str.extract('gene=([\w-]*);')
# REFMAP substring extracted from qry_id_list column
REFMAPgeneIDs = REFMAP['qry_id_list'].str.extract('([\w-]*)|')

# index GFF with gene IDs
GFF.index = GFFgeneIDs[0]

# make set of REFMAP geneIDs
REFMAPgeneIDset = set(REFMAPgeneIDs[0])

# loop through geneIDs in GFF, skipping nan
for GENE in set(GFF.index):
	if pd.notna(GENE):
		
		# skip features if geneID found in REFMAP
		if GENE not in REFMAPgeneIDset:
			
			# append all features from geneID to OUTGFF
			GFF.loc[GENE].to_csv(OUTGFF, sep='\t',
								quoting=csv.QUOTE_NONE,
								index=False, header=False, 
								mode='a')
	















