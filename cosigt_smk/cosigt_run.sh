#!/bin/bash
set -x

#conda environemnt - load
eval "$(conda shell.bash hook)"
module load singularity/3.8.5
conda activate snakemakeenv_latest
#run
bindings=$(cat singularity_bind_paths.csv)
stringb=$(echo "-B $bindings")
snakemake evaluation --use-singularity --singularity-args "$stringb" --rerun-triggers mtime --profile config/slurm
