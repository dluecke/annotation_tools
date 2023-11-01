#!/usr/bin/env python3

# duplicate_regions_from_BLAST.py takes BLAST files for self-search and sister-search
# returns list of regions with better hits in same annotation as in sister assembly
# written as part of pipeline which searches multifasta of gene regions against same

import pandas as pd

INFILE_WITHIN = "AllVsAll-blast_results-min1kb.tsv"
INFILE_BETWEEN = "AllVsrhop-blast_results-min1kb.tsv"

# parameters for reading BLAST outfmt 6 tsv files
outfmt6_kept_columns_pos = [0, 1, 2, 3, 11]
outfmt6_kept_columns_names = ['Query', 'Subject', 'PctID', 'Length', 'BitScore']

# read BLAST output from within assembly search
df_within = pd.read_csv(INFILE_WITHIN, delimiter='\t',
	usecols=outfmt6_kept_columns_pos, names=outfmt6_kept_columns_names)
# remove match-to-self lines
df_within = df_within[df_within['Query'] != df_within['Subject']]
# label as not between species
df_within['Between'] = False

# read BLAST output from search to sister assembly
df_between = pd.read_csv(INFILE_BETWEEN, delimiter='\t',
	usecols=outfmt6_kept_columns_pos, names=outfmt6_kept_columns_names)
# label as between species
df_between['Between'] = True

# build full df of best hit (by BitScore) matches with 2-level index
df_matches = pd.concat(
	[df_within, df_between], ignore_index=True
	).sort_values(by='BitScore', ascending=False)
# sort by BitScore and only keep top for a given Query-Subject pair
df_matches = df_matches.drop_duplicates(subset=['Query', 'Subject'])
# set index as unique Query-Subject pairs
df_matches.set_index(['Query', 'Subject'], inplace=True)
# ensure BitScore sorting within index structure
df_matches = df_matches.sort_index().sort_values(by='BitScore', ascending=False)

# call duplicates as all within-assembly matches with higher BitScore than any between match
df_duplicates = df_matches.groupby(level='Query', group_keys=False).apply(
	lambda query: query.loc[query['Between'].cumsum() < 1]
)
# drop index
df_duplicates = df_duplicates.reset_index()

# get reciprocally called duplicates
df_duplicatesRecip = df_duplicates[df_duplicates['Query'].isin(df_duplicates['Subject'])]

# extract gene names, scaffold positions from reciprocal duplicate matches
df_regionsQ = df_duplicatesRecip['Query'].str.extract(r'(.*)::(Scaffold_.*):(.*)-(.*)\(.*')
df_regionsQ.columns = ['Q-GeneID', 'Q-Scaffold', 'Q-Start', 'Q-Stop']
df_regionsS = df_duplicatesRecip['Subject'].str.extract(r'(.*)::(Scaffold_.*):(.*)-(.*)\(.*')
df_regionsS.columns = ['S-GeneID', 'S-Scaffold', 'S-Start', 'S-Stop']
# sort with same method as geneIDs were assigned
df_regions = pd.concat([df_regionsQ, df_regionsS], axis=1).sort_values(by=['Q-Scaffold', 'Q-Start'])

# get list of genes with reciprocal duplicate matches
df_genes = pd.DataFrame()
df_genes['GeneID'] = df_regions['Q-GeneID'].append(df_regions['S-GeneID']).str.split(',')
df_genes = df_genes.explode('GeneID').drop_duplicates().reset_index(drop=True)
df_genes['GeneNum'] = df_genes['GeneID'].str.extract(r'ID=gene(.*)').astype(int)
df_genes = df_genes.sort_values(by='GeneNum').reset_index(drop=True)

# want all scaffold pairs that contain putative duplicates, without reverse pairs
set_DupScaffTuples = set()
# sorted list of the tuple set
ls_DupScaffTuples_wRPs = sorted(set(zip(df_regions['Q-Scaffold'], df_regions['S-Scaffold'])))
# add tuples to scaffold pair set, skipping second member of reverse pairs
for ScaffPair in ls_DupScaffTuples_wRPs:
	if (ScaffPair[1], ScaffPair[0]) not in set_DupScaffTuples:
		set_DupScaffTuples.add(ScaffPair)
# ordered list of unique putative duplicated scaffold pairs
ls_DupScaffs = sorted(set_DupScaffTuples)

# make df for coordinates of IDed duplication for each pair
df_DupCoords = pd.DataFrame()
for SP in ls_DupScaffs:
	SPDCs = df_regions[ (df_regions['Q-Scaffold']==SP[0]) & (df_regions['S-Scaffold']==SP[1]) ]







