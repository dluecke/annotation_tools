# annotation_tools
Assorted tools for working with genome annotation and assembly qc

## Used in *D. prolongata* genome report Luecke *et al* 2024

### Combining GFF files from multiple sources
The `merge_GFFs/` directory has scripts used to extract the novel features in a new GFF not present in a reference file.

### Identifying and removing wholly duplicate scaffolds
The `remove_duplicate_scaffolds` directory has scripts to parse BUSCO output and sort scaffolds by proportion of "Duplicate" BUSCO benchmarks, then to split verified duplicate scaffolds from assembly and annotation.

### Identifying candidate duplicate genes
The `id_duplicate_genes/` directory has scripts to extract all annotated gene regions from assembly/annotation into multi-fasta for reciprocal BLAST search, then to parse the resulting BLAST output to flag gene regions wtih higher similarity to other regions of same assembly than to sister species reference assembly.

### Whole genome alignment and visualization with mummer
The `alignment_and_visualization/` directory has scripts used in conjunction with the mummer aligner, to extract priority scaffolds before running alignment and to increase flexibility of plotting for alignment output.


## Other tools

### Scoring annotation for gene family coverage
The `score_annotations_by_BLAST/` directory has scripts to compare assembly regions with BLAST hits to regions with identified gene models.
