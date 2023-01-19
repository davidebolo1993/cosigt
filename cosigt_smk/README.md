# Snakemake pipeline

## prepare config and samples.tsv, organize resources

```bash
./cosigt_prepare.sh <path/to/cram_directory> <path/to/reference.fa> <path/to/graph.gfa> <region> <blacklist.txt> #chr1:103456064-103863972 for AMY, for instance
```

## run the snakemake pipeline
Uncommenting --profile will run on slurm cluster

```bash
#assume snakemake and singularity in path for cosigt, R with ggplot2/data.table/rjson in path for evaluate
./cosigt_run.sh
```

