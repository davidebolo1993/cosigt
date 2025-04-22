#!/bin/bash
SINGULARITY_TMPDIR=/scratch/davide.bolognini snakemake --profile config/slurm --singularity-args "-B /scratch/davide.bolognini,/group/soranzo/davide.bolognini/working/dev/cosigt_paper/real_data -e" cosigt
