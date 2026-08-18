# Cosigt snakemake pipeline

Docs for the pipeline are available [online](https://davidebolo1993.github.io/cosigtdoc/).

## Quickstart

There are three commands. From the repository root:

```bash
make init     # create config/ from the examples and write .cosigt.mk
make check    # validate the environment and the configuration
make run      # run the pipeline
```

Edit the files created in `cosigt_smk/config/` between `init` and `check`.

### make check

`check` validates everything needed for the settings you give it, and fails with
a specific message rather than part-way through a run. It checks:

- Snakemake is on PATH.
- The executor plugin required by `PROFILE` is installed (`slurm`, `lsf`,
  `cluster-generic`); for `local` it reports the cores it will use.
- The runtime required by `SOFTWARE` is available: `apptainer` (accepting either
  `apptainer` or `singularity`), `conda` (including the 24.7.1 minimum Snakemake
  needs), or, for `none`, that every tool the current configuration will invoke
  is on PATH. That tool list is derived from the config, so it accounts for
  `read_mode`, `vcf`, `gtf` and the visualisation switches.
- The config file, sample table, region BED, assembly or allele tables, and all
  referenced files and indexes, writing the generated region metadata and the
  normalised flagger blacklist as a side effect.
- Finally, it composes the Apptainer flags this run needs and writes them to
  `cosigt_smk/.cosigt/apptainer.args`.

Those flags are the bind mounts covering every configured input and output
location, collapsed to the shortest set of parent directories, plus `-e`
(`--cleanenv`), which pggb requires. `make run` picks the file up automatically
. Set `apptainer_extra` in
`config.yaml` to append site-specific flags, or `apptainer_cleanenv: false` to
drop `-e`.

### make run

`run` just runs the selected target. It re-reads the flags written by `check`,
so run `check` first, or after changing any input path.

### Settings

All three commands take the same variables, on the command line or persisted in
`.cosigt.mk` by `make init`:

| variable | default | meaning |
| --- | --- | --- |
| `PROFILE` | `local` | `local`, `slurm`, `lsf`, `cluster-generic`, a path, or `none` |
| `SOFTWARE` | `apptainer` | `apptainer`, `conda`, or `none` |
| `TARGET` | `cosigt` | Snakemake target: `cosigt`, `refine`, `benchmark` |
| `CORES` | all detected | passed to `--cores` |
| `SMK_ARGS` | empty | extra Snakemake arguments, e.g. `-n` for a dry run |

So a cluster run is `make run PROFILE=slurm`, a dry run is
`make run SMK_ARGS=-n`, and a site-specific submit command is
`make run PROFILE=cluster-generic SMK_ARGS='--cluster-generic-submit-cmd "qsub"'`.

### Benchmarking

`make run TARGET=benchmark` compares the haplotypes cosigt predicts against each
sample's true haplotypes. Predicted and true sequences are aligned with edlib and
reported as QV and error rate. `benchmark_mode` selects the experiment:

**`leave_zero_out`** (default). The genotyped samples' own haplotypes are part of
the graph, so an exact hit is possible and QV measures whether cosigt finds it.
This is an upper bound on performance.

**`leave_all_out`**. The graph cosigt genotypes against contains no haplotype from
any genotyped sample, so it must reconstruct each sample from other people's
haplotypes. Build that graph by pointing `assemblies` at a panel that excludes the
genotyped samples. The truth then has to come from elsewhere, so this mode also
needs `truth_graphs`: a TSV mapping each region to a pre-built odgi graph that
*does* contain the samples' own haplotypes.

```
region                graph
chr1_100000_120000    /path/to/chr1_100000_120000.truth.og
chr2_200000_230000    /path/to/chr2_200000_230000.truth.og
```

Both predicted and true sequences are read from that graph, so they share one
orientation frame. Because a perfect answer is impossible by construction, each
true haplotype is also compared against every haplotype the genotyping graph did
offer:

```
QV_sum_best = best(hap_1_true vs panel) + best(hap_2_true vs panel)
QVfrac     = QV_sum_pred / QV_sum_best
```

`QVfrac` near 1 means cosigt chose about the best pair available to it, so a low
absolute QV reflects an incomplete panel rather than a genotyping error.

`make check` verifies that the table is present, covers every region, and that
each graph exists and is an `.og` (convert a GFA with
`odgi build -g in.gfa -o out.og`). Note that this scan is the expensive part of
the mode: it runs one edlib alignment per (true haplotype, panel haplotype) pair,
so raise `benchmark.runtime` for large panels.

The per-region tables are collected into `benchmark/benchmark.qv.tsv`:

| column | meaning |
| --- | --- |
| `sample`, `region`, `gene_name` | identifiers; `gene_name` is column 4 of the regions BED |
| `hap_1_pred`, `hap_2_pred` | predicted haplotypes |
| `cluster_1_pred`, `cluster_2_pred` | their clusters, as assigned by `cluster.r` |
| `hap_1_true`, `hap_2_true` | the sample's own haplotypes, in the matched order |
| `QV_1_pred`, `QV_2_pred`, `QV_sum_pred` | phred-scaled quality per haplotype, and their sum |
| `error_rate_1_pred`, `error_rate_2_pred` | edit distance over alignment length |
| `hap_1_best`, `hap_2_best` | best-matching haplotype available in the genotyping graph (`leave_all_out` only) |
| `QV_1_best`, `QV_2_best`, `QV_sum_best` | the QV those would have achieved |
| `QVfrac` | `QV_sum_pred / QV_sum_best` |

Alongside the table, `benchmark/benchmark.qv.png` summarises the result per gene
as a stacked bar, genes ordered by the share in the best band, with the number of
samples above each bar. What is banded depends on the mode:

- `leave_zero_out` bands `QV_pred` into four quality bands (`<=17`, `17-23`,
  `23-33`, `>33`), ordered by the share reaching `>33`.
- `leave_all_out` bins `QVfrac` into quintiles, ordered by the share in
  `Q5: 0.8-1.0`, since there the achieved QV alone cannot separate genotyping
  error from a panel that held nothing closer.

`benchmark.max_bars_per_row` (default 30) splits the plot into rows.

Both assignments of the two predictions to the two truths are scored and the
better one is reported, so haplotype order does not matter. Rows are emitted with
`NA` values rather than dropped when a sample's own haplotypes are not in the
graph exactly twice, and when a region's ploidy is not 2, since the comparison is
defined for two predicted against two true haplotypes. Note that `compute_qv`
floors the edit distance at 0.5, so an exact match gives a length-dependent QV
ceiling rather than infinity — about 36 for a 2 kb haplotype, 53 for 100 kb.

This target is container-only: `compute_qv` has no conda package, so run
benchmarking with `SOFTWARE=apptainer`.

## Job grouping on clusters

The per-sample, per-region part of the workflow is a linear chain of about seven
rules (read extraction, optional k-mer rescue, realignment, graph injection,
coverage, genotyping). Submitted individually that is `7 x samples x regions`
cluster jobs, which dominates the total: the whole workflow is roughly
`7*S*R + 15*R + S` jobs, so 500 samples over 100 regions is around 350k
submissions.

Those rules are therefore assigned to a Snakemake group called `genotype`.
Snakemake submits each connected component of a group as a single cluster job
and runs the rules inside it in order, which cuts those submissions by 7x with
no change to the DAG, to rule granularity, or to resumability. Each rule still
uses its own container or conda environment inside the group job.

`--group-components genotype=N` packs N independent (sample, region) chains into
one submission, giving a dial between full parallelism and a few large jobs. It
is set in the cluster profiles and defaults to 1:

```bash
make run PROFILE=slurm SMK_ARGS='--group-components genotype=8'
```

Group resources are derived automatically, so no group-level resources need to
be declared. Resources of rules that run in series are combined by taking the
maximum, except `runtime`, which is summed; rules that can run in parallel have
their resources summed and their `runtime` maximised. With the default config a
single `genotype` group job requests about 41 minutes, 22 GB and 10 cores, the
memory and cores being driven by `kfilt` and `samtools_fasta_mapped` sharing the
first layer. Raising `--group-components` multiplies the packed chains, so
increase walltime accordingly if jobs start hitting limits.

Set `allele_source` in `config.yaml` to `assemblies` for the standard workflow
or `custom` for user-provided per-region allele FASTAs. Set `read_mode` to
`short`, `ancient`, or `long:<preset>` to choose the read realignment strategy.
Supported long-read minimap2 presets are `map-pb`, `map-hifi`, `map-ont`,
`map-iclr`, and `lr:hq`; for example, use `read_mode: long:map-ont` for ONT
reads.

Unmapped-read rescue (meryl + kfilt) runs in `short` mode only. There, reads
that failed to map to the reference are recovered by matching region-specific
31-mers and realigned along with the mapped ones. It is deliberately skipped for
`long:<preset>` and `ancient`, where per-read error rates make exact 31-mer
matching an unreliable filter. In those modes only reads mapped within the
region are realigned, and the `meryl` and `kfilt` resource blocks in
`config.yaml` are unused.

The workflow keeps reusable graph, allele FASTA, aligner index, region metadata,
sample unmapped-read FASTA, and final genotype outputs. Large per-sample,
per-region transport files such as mapped FASTA slices, filtered FASTA,
realigned CRAMs, GAF, and gafpack coverage are temporary. This keeps reruns
incremental: adding samples reuses existing region graphs and indexes, while
adding regions reuses chromosome-level assembly mappings, reference k-mer data,
and sample-level unmapped reads.