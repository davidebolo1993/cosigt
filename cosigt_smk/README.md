# Cosigt snakemake pipeline

How to run the pipeline. Everything else — configuration keys, benchmarking, cluster tuning, worked examples — is in the
[online docs](https://davidebolo1993.github.io/cosigtdoc/).

## Quickstart

There are three commands. From the repository root:

```bash
make init     # create config/ from the examples and write .cosigt.mk
make check    # validate the environment and the configuration
make run      # run the pipeline
```

Edit the files created in `cosigt_smk/config/` between `init` and `check`.

### make init

`init` creates `cosigt_smk/config/` and copies each shipped example to the name the pipeline actually reads:

| created | what goes in it |
| --- | --- |
| `config.yaml` | every pipeline setting |
| `samples.tsv` | one row per sample to genotype |
| `regions.bed` | regions of interest, with per-region ploidy |
| `assemblies.tsv` | assemblies per chromosome, for `allele_source: assemblies` |
| `alleles.tsv` | curated panel, for `allele_source: alleles` |

Only one of `assemblies.tsv` and `alleles.tsv` is used, depending on `allele_source`. `truth_graphs.tsv.example` is not copied, because it is only needed for `benchmark_mode: leave_all_out`; copy it yourself if you need it.

It also writes `.cosigt.mk` with the settings below, so later commands can be run without repeating them. Existing files are never overwritten, so `init` is safe to re-run and will not touch edits you have already made.

### make check

`check` validates everything needed for the settings you give it, and fails with a specific message rather than part-way through a run. It checks:

- Snakemake is on PATH.
- The executor plugin required by `PROFILE` is installed (`slurm`, `lsf`, `cluster-generic`); for `local` it reports the cores it will use.
- The runtime required by `SOFTWARE` is available: `apptainer` (accepting either `apptainer` or `singularity`), `conda` (including the 24.7.1 minimum Snakemake needs), or, for `none`, that every tool the current configuration will invoke is on PATH. That tool list is derived from the config, so it accounts for `read_mode`, `vcf`, `gtf` and the visualisation switches.
- The config file, sample table, region BED, assembly or allele tables, and all referenced files and indexes, writing the generated region metadata and the normalised flagger blacklist as a side effect.
- Finally, it composes the Apptainer flags this run needs, writes them to `cosigt_smk/.cosigt/apptainer.args`, and persists the validated settings to `.cosigt.mk`.

Those flags are the bind mounts covering every configured input and output location, collapsed to the shortest set of parent directories, plus `-e` (`--cleanenv`), which pggb requires. `make run` picks the file up automatically.
Set `apptainer_extra` in `config.yaml` to append site-specific flags, or `apptainer_cleanenv: false` to drop `-e`.

### make run

`run` just runs the selected target. It re-reads the flags written by `check`, so run `check` first, or after changing any input path.

### Settings

All three commands take the same variables, on the command line or persisted in `.cosigt.mk` by `make init` and `make check`:

| variable | default | meaning |
| --- | --- | --- |
| `PROFILE` | `local` | `local`, `slurm`, `lsf`, `cluster-generic`, a path, or `none` |
| `SOFTWARE` | `apptainer` | `apptainer`, `conda`, or `none` |
| `TARGET` | `cosigt` | Snakemake target: `cosigt`, `graph`, `refine`, `benchmark` |
| `CORES` | all detected | passed to `--cores`; on a cluster, size it from a compute node |
| `SMK_ARGS` | empty | extra Snakemake arguments, e.g. `-n` for a dry run |

`TARGET=graph` builds the per-region pangenome graphs, their clustering and the
optional plots, and stops there: nothing that reads a sample's alignment runs.
Useful to get the expensive graph construction done once, or in advance of
having the reads. A later `TARGET=cosigt` reuses everything it produced.

So a cluster run is `make run PROFILE=slurm`, a dry run is
`make run SMK_ARGS=-n`, and a site-specific submit command is
`make run PROFILE=cluster-generic SMK_ARGS='--cluster-generic-submit-cmd "qsub"'`.
