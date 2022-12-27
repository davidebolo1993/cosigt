# Snakemake pipeline

## prepare config and samples.tsv, organize resources

```bash
./cosigt_prepare.sh <path/to/cram_directory> <path/to/reference.fa> <path/to/graph.gfa> <region> #"chr1:103456064-103863972 for AMY, for instance
```

## run the snakemake pipeline
Uncommenting --profile will run on slurm cluster, but first ajust the slurm config

```bash
#edit config/slurm/config.yaml first 
./cosigt_run
```

