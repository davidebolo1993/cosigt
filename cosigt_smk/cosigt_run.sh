#!/bin/bash
set -x

#conda environemnt - load
eval "$(conda shell.bash hook)"
conda activate /global/home/users/davidebolognini/micromamba/envs/snakemake_latest
#run
bindings=$(cat singularity_bind_paths.csv)
stringb=$(echo "-B $bindings")
snakemake --unlock
snakemake --cores 5 --use-singularity --singularity-args "$stringb" cosigt #--profile config/slurm
snakemake evaluate --cores 1 --use-conda #--profile config/slurm
