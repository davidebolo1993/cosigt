#!/bin/bash
set -x

#conda environemnt - load
eval "$(conda shell.bash hook)"
conda activate /global/home/users/davidebolognini/micromamba/envs/snakemake_latest
#run
bindings=$(cat singularity_bind_paths.csv)
stringb=$(echo "-B $bindings")
snakemake --unlock
snakemake --cores 10 --use-singularity --singularity-args "$stringb" cosigt #--profile config/slurm
module load r/3.6.3
module load r-packages/default
snakemake --cores 1 evaluate #--use-conda to try to install the r env
