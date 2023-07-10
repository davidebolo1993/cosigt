# Snakemake pipeline

```bash
#blacklist is optional
python workflow/scripts/prepare.py  <path/to/cram_directory> <path/to/reference.fa> <path/to/graph.gfa> <region> <blacklist.txt>
#assume snakemake and singularity in path
#run evaluation
bindings=$(cat singularity_bind_paths.csv)
stringb=$(echo "-B $bindings")
snakemake all --use-singularity --singularity-args "$stringb" --rerun-triggers mtime --profile config/slurm
```


