# Cosigt snakemake pipeline

Docs for the pipeline are available [online](https://davidebolo1993.github.io/cosigtdoc/).

## Unified workflow quickstart

From the repository root:

```bash
make init
```

Edit the files created in `cosigt_smk/config/`, then validate them without
running compute jobs:

```bash
make check
```

`make check` only builds the small validation target. It checks the config files,
sample table, region BED, indexes, required input files, and generated metadata
paths. Use it as the fast "is my setup sane?" command.

Run a dry-run of the full genotyping target:

```bash
make dryrun
```

`make dryrun` builds the complete DAG for the selected target, `cosigt` by
default, but still does not run jobs. Use it to see which real pipeline steps
would run and whether all wildcards and dependencies resolve.

Run locally:

```bash
make run
```

The Makefile uses Snakemake's conda environments by default. On systems with
Apptainer/Singularity available, add `SOFTWARE=apptainer` to any run command.
With current Snakemake versions, `SOFTWARE=conda` requires conda 24.7.1 or
newer in the environment that runs Snakemake. The Makefile checks this before
launching Snakemake and prints a short update hint if conda is too old.

Run with the SLURM executor plugin:

```bash
make run-slurm
```

Run with the LSF executor plugin:

```bash
make run-lsf
```

The Makefile checks for the relevant Snakemake executor plugin before cluster
runs. If one is missing, install just that plugin with the command it prints, or
install the common cluster plugins together:

```bash
make install-cluster-plugins
```

Set `allele_source` in `config.yaml` to `assemblies` for the standard workflow
or `custom` for user-provided per-region allele FASTAs. Set `read_mode` to
`short`, `ancient`, `ont`, `pacbio_hifi`, or `pacbio_clr` to choose the read
realignment strategy.

The workflow keeps reusable graph, allele FASTA, aligner index, region metadata,
sample unmapped-read FASTA, and final genotype outputs. Large per-sample,
per-region transport files such as mapped FASTA slices, filtered FASTA,
realigned CRAMs, GAF, and gafpack coverage are temporary. This keeps reruns
incremental: adding samples reuses existing region graphs and indexes, while
adding regions reuses chromosome-level assembly mappings, reference k-mer data,
and sample-level unmapped reads.
