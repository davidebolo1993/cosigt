#!/usr/bin/env bash
# Snakemake runs this automatically after creating the r.yaml conda environment.
# SVbyEye is not packaged on conda-forge/bioconda, so it is installed from
# GitHub exactly as the davidebolo1993/renv:4.3.3 container does. Without it,
# the plot_ava rule (workflow/scripts/plotava.r) fails under
# --software-deployment-method conda.
set -euo pipefail

Rscript -e "devtools::install_github('daewoooo/SVbyEye', branch='master', upgrade='never')"
