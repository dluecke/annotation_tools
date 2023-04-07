#!/usr/bin/env python3

# write_novel_gff_features.py extracts all features marked as novel by gffcompare reference comparison
# takes ($1) NEW.gff annotation and $(2) gffcmp.NEW.gff.tmap 
# may need to strip NEW.gff of non-feature lines for MAKER (written using Liftoff output)
# takes qry_gene_id from features with class_code i, y, p, or u from tmap
#    and extracts corresponding lines from NEW.gff
# outputs NEW.gff-unique.gff (these lines can be added to reference gff without introducing redundancy)
# use gffcompare -r REFERENCE.gff NEW.gff to get tmap file


import sys
import pandas as pd
import csv

GFF = pd.read_csv(sys.argv[1], sep='\t', header=None)

TMAP = pd.read_csv(sys.argv[2], sep='\t', index_col=2)

OUTGFF = sys.argv[1] + "-unique.gff"

# get i, y, p, and u transcripts from TMAP
UniqueKeys = ['i', 'y', 'p', 'u']

TMAPunique = pd.DataFrame()
for KEY in UniqueKeys:
	TMAPunique = TMAPunique.append(TMAP.loc[KEY])

# set of unique geneIDs 
TMAPuniqueGenes = set(TMAPunique['qry_gene_id'])


# index GFF by geneID
# GFF gene=...; substring extracted from 9th column
GFFgeneIDs = GFF[8].str.extract('gene=([\w-]*);')
# some Liftoff uses gene_name= instead of gene=, use this if gene= gives NA
GFFgeneIDs['gene_name'] = GFF[8].str.extract('gene_name=([\w-]*);')
GFFgeneIDs[0].fillna(GFFgeneIDs.gene_name, inplace = True)
# index GFF with gene IDs
GFF.index = GFFgeneIDs[0]

# print all GFF features from unique gene set 
for GENE in TMAPuniqueGenes:

	# use double brackets to ensure dataframe structure (otherwise single rows become series)
	GFF.loc[[GENE]].to_csv(OUTGFF, sep='\t',
						quoting=csv.QUOTE_NONE,
						index=False, header=False, 
						mode='a')


