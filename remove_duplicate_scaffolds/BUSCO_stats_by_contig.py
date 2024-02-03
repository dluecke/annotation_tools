#!/usr/bin/env python3

# BUSCO_stats_by_contig.py reads the BUSCO program output file "full_table.tsv"
# returns table of contig breakdown and list of fully duplicated contigs for removal

# sorts by contig/scaffold ("Sequence") in the original assembly and reports for each:
# 	N_BUSCOs: total BUSCO benchmark genes (BUSCOs)
# 	n_Complete: number "Complete" BUSCOs
#	p_Complete: percent BUSCOs which are "Complete"
#	n_Duplicated: number "Duplicated" BUSCOs
#	p_Duplicated: percent BUSCOs which are "Duplicated"
#	n_Fragmented: number "Fragmented" BUSCOs
#	p_Fragmented: percent BUSCOs which are "Fragmented"
#	n_InternalDup: number of duplicates from within the contig eg tandem duplicates
#	"Seqs w Duplicate Pair": All other contigs with copies of Duplicated BUSCOs
#	"Duplicated BUSCO IDs": IDs for all Duplicated BUSCOs

# Checks all contigs with only duplicated BUSCOs for removal and outputs 
# (if a BUSCO is only present on multiple all-duplicates contigs need to keep one)
# After check reports number of contigs and BUSCO duplicates which can be removed

# writes .tsv tables named as modifications of input tsv INFILE.tsv
#	INFILE_SeqStats.tsv: stats for all contigs
#	INFILE_SeqStats-Keep.tsv: subset of contigs with unique BUSCOs
#	INFILE_SeqStats-Remove.tsv: subset of contigs with only duplicated BUSCOs
#	INFILE_RemoveList.tsv: just the suggested contigs to Remove and the corresponding Keep contigs

import sys
import pandas as pd

# Read in BUSCO table as CL argument 1, skipping 2 comment lines at top
INFILEext = sys.argv[1]
BUSCO_full_table = pd.read_csv(INFILEext, delimiter='\t', skiprows=2)

# Prep output filenames
INFILE = INFILEext[:INFILEext.rfind('.')]
OUTFILE_ALL = INFILE + "_SeqStats.tsv"
OUTFILE_KEEP = INFILE + "_SeqStats-Keep.tsv"
OUTFILE_REMOVE = INFILE + "_SeqStats-Remove.tsv"
OUTFILE_REMOVELIST = INFILE + "_RemoveList.tsv"

# Extract only rows for BUSCOs found in assembly (Complete, Duplicated, or Fragmented)
df_FoundBUSCOs = BUSCO_full_table[BUSCO_full_table['Status'] != 'Missing']
df_FoundBUSCOs = df_FoundBUSCOs.rename(columns={df_FoundBUSCOs.columns[0]: 'Busco ID'})

# Function to find all BUSCOs for given scaffold and calculate statistics
def ContigStats(ROW):
	# extract only rows for given Sequence
	DF_CONTIG = df_FoundBUSCOs[df_FoundBUSCOs['Sequence'] == ROW['Sequence']]
	ROW['N_BUSCOs'] = DF_CONTIG.shape[0]
	ROW['n_Complete'] = sum(DF_CONTIG['Status'] == 'Complete')
	ROW['p_Complete'] = ROW['n_Complete'] / ROW['N_BUSCOs']
	ROW['n_Duplicated'] = sum(DF_CONTIG['Status'] == 'Duplicated')
	ROW['p_Duplicated'] = ROW['n_Duplicated'] / ROW['N_BUSCOs']
	ROW['n_Fragmented'] = sum(DF_CONTIG['Status'] == 'Fragmented')
	ROW['p_Fragmented'] = ROW['n_Fragmented'] / ROW['N_BUSCOs']
	ROW['n_InternalDup'] = ROW['N_BUSCOs'] - len(set(DF_CONTIG['Busco ID']))
	# extract all rows with shared BUSCOs to report scaffold pairs
	DF_ALLDUPS = df_FoundBUSCOs[df_FoundBUSCOs['Busco ID'].isin(DF_CONTIG['Busco ID'])]
	DUP_CONTIGS = set(DF_ALLDUPS['Sequence'])
	DUP_CONTIGS.discard(ROW['Sequence'])
	ROW['Seqs w Duplicate Pair'] = DUP_CONTIGS
	# set of all Duplicated BUSCOs in Sequence
	ROW['Duplicated BUSCO IDs'] = set(DF_CONTIG[DF_CONTIG['Status'] == 'Duplicated']['Busco ID'])
	return(ROW)

# Make and fill dataframe for statistics
df_StatsByScaffold = pd.DataFrame({'Sequence': sorted(set(df_FoundBUSCOs['Sequence']), key = lambda x: (len(x), x))})	
df_StatsByScaffold = df_StatsByScaffold.apply(ContigStats, axis=1)
df_StatsByScaffold = df_StatsByScaffold.set_index('Sequence')

# Split stats dataframe into RemoveCandidates if all BUSCOs are Duplicated, and Keepers otherwise
df_RemoveCandidates = df_StatsByScaffold[df_StatsByScaffold['p_Duplicated'] == 1].sort_values('N_BUSCOs')
df_Keepers = df_StatsByScaffold[df_StatsByScaffold['p_Duplicated'] < 1]

# Initialize empty dataframe for Remove set
df_Remove = pd.DataFrame(index=pd.Index([], name='Sequence'), columns=df_StatsByScaffold.columns)

# Iterate through RemoveCandidates 
for INDEX, ROW in df_RemoveCandidates.iterrows():
	ROW_DUP_IDS = df_RemoveCandidates.loc[INDEX]['Duplicated BUSCO IDs']
	KEPT_DUP_IDS = set().union(*df_Keepers['Duplicated BUSCO IDs'].values)
	# add to Keep set if ROW has BUSCO IDs not already in Keepers
	if len(ROW_DUP_IDS & KEPT_DUP_IDS) < len(ROW_DUP_IDS):
		df_Keepers = df_Keepers.append(ROW)
	# otherwise add to Remove 
	else:
		df_Remove = df_Remove.append(ROW)
		
# print removal summary results to screen
N_REMOVE_CONTIGS = df_Remove.shape[0]
N_REMOVE_BUSCOS = len(set().union(*df_Remove['Duplicated BUSCO IDs'].values))
print()
print(N_REMOVE_CONTIGS, "sequences contain only duplicated BUSCOs")
print("Removing these will remove duplicate copies from", N_REMOVE_BUSCOS, "BUSCOs")
print("-----------")

# Initialize 2 column dataframe for Remove scaffolds and corresponding Keep scaffolds
df_RemoveList = pd.DataFrame(index=pd.Index(sorted(set(df_Remove.index), key = lambda x: (len(x), x)), name='Remove'), columns=['Keep'] )
for i in df_RemoveList.index:
	# list of scaffolds with BUSCO pairs not being removed
	KEEP = ', '.join(df_Remove.loc[i]['Seqs w Duplicate Pair'] - set(df_Remove.index))
	df_RemoveList.loc[i]['Keep'] = KEEP

# write output to files
df_StatsByScaffold.to_csv(OUTFILE_ALL, sep='\t')
df_Keepers.to_csv(OUTFILE_KEEP, sep='\t')
df_Remove.to_csv(OUTFILE_REMOVE, sep='\t')
df_RemoveList.to_csv(OUTFILE_REMOVELIST, sep='\t')

# print outfile descriptions to screen
print("Files with statistics for BUSCOs by sequence:", OUTFILE_ALL, OUTFILE_KEEP, OUTFILE_REMOVE, sep='\n')
print("Containing all sequences, sequences with unique BUSCOs, and sequences with only duplicated BUSCOs respectively.")
print("\nFile with list of scaffolds suggested for removal and scaffolds containing corresponding unremoved BUSCO copies:")
print(OUTFILE_REMOVELIST, "\n")




