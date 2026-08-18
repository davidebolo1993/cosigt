import csv
import os
import shlex
from collections import OrderedDict
from shutil import copyfile, which

from snakemake.exceptions import WorkflowError
from snakemake.utils import validate


WORKDIR = os.getcwd()
SCHEMA_DIR = os.path.join(WORKDIR, "workflow", "schemas")
CONFIG_SCHEMA = os.path.join(SCHEMA_DIR, "config.schema.yaml")
SAMPLES_SCHEMA = os.path.join(SCHEMA_DIR, "samples.schema.yaml")
ASSEMBLIES_SCHEMA = os.path.join(SCHEMA_DIR, "assemblies.schema.yaml")
ALLELES_SCHEMA = os.path.join(SCHEMA_DIR, "alleles.schema.yaml")
TRUTH_GRAPHS_SCHEMA = os.path.join(SCHEMA_DIR, "truth_graphs.schema.yaml")


RESOURCE_DEFAULTS = {
    # max_occ is bwa-mem2's -c. It aborts on an internal assert when too many
    # seeds share a position, which an allele panel provokes by construction;
    # lowering -c is the documented workaround, so each retry halves it down to
    # min_max_occ. See bwa-mem2 issue #269.
    "bwa-mem2": {
        "threads": 5, "mem_mb": 10000, "runtime": 4,
        "max_occ": 500, "min_max_occ": 50, "retries": 3,
    },
    "bwa": {"threads": 5, "mem_mb": 10000, "runtime": 4},
    "minimap2": {
        "ava": {"threads": 4, "mem_mb": 5000, "runtime": 2},
        "avo": {"threads": 6, "mem_mb": 40000, "runtime": 50},
        "reads": {"threads": 8, "mem_mb": 16000, "runtime": 8},
    },
    "samtools": {
        "fasta_mapped": {"threads": 2, "mem_mb": 2000, "runtime": 2},
        "fasta_unmapped": {"threads": 8, "mem_mb": 8000, "runtime": 8},
    },
    "pggb": {
        "threads": 24,
        "mem_mb": 20000,
        "runtime": 40,
        "tmpdir": "/tmp",
        "params": "-c 2 -k 101",
    },
    # pangene_prepare runs one miniprot per allele contig, in parallel.
    "pangene": {"threads": 8},
    "impg": {"threads": 4},
    "panplexity": {"threads": 4},
    "meryl": {"threads": 10, "mem_mb": 25000, "runtime": 20},
    "kfilt": {"threads": 8, "mem_mb": 20000, "runtime": 20},
    # benchmark_qv runs one edlib alignment per (true haplotype, panel
    # haplotype) pair, so it scales with panel size x samples and is worth
    # parallelising.
    "benchmark": {"threads": 8, "mem_mb": 10000, "runtime": 60, "max_bars_per_row": 30},
    "default": {
        "small": {"mem_mb": 500, "runtime": 2},
        "mid": {"mem_mb": 2000, "runtime": 5},
        "high": {"mem_mb": 10000, "runtime": 10},
    },
}


BOOLEAN_DEFAULTS = {
    "wally_viz": False,
    "svbyeye_viz": False,
    "odgi_viz": True,
    "pangene_viz": True,
    "sv_calling": False,
    "vcf": False,
}


# Tuning knobs for workflow/scripts/cluster.r. Keys mirror the script's
# --kebab-case options and are rendered back into that form by cluster_args().
CLUSTER_DEFAULTS = {
    "similarity_threshold": "automatic",
    "levels": 1,
    "eps_selection": "quality",
    "eps_min": 0,
    "eps_max": 0.30,
    "eps_step": 0.01,
    "score_use_dbcv": True,
    "dbcv_dim": 2,
    "low_diversity_mpd_norm": 0.05,
    "giant_cluster_fraction": 0.85,
    "small_cluster_size": 2,
    "ambiguous_eps_ratio": 1.25,
    "plot_clusters": False,
    "plot_label_max_cluster_size": 1,
    "plot_width": 8,
    "plot_height": 6,
}

CLUSTER_BOOLEAN_KEYS = ("score_use_dbcv", "plot_clusters")


def _fail(message):
    raise WorkflowError(message)


def _deep_defaults(target, defaults):
    for key, value in defaults.items():
        if isinstance(value, dict):
            target.setdefault(key, {})
            if not isinstance(target[key], dict):
                _fail(f"Config key '{key}' must be a mapping.")
            _deep_defaults(target[key], value)
        else:
            target.setdefault(key, value)


def _migrate_time_to_runtime(node):
    if isinstance(node, dict):
        if "time" in node and "runtime" not in node:
            node["runtime"] = node["time"]
        for value in node.values():
            _migrate_time_to_runtime(value)


def _normalize_bool(value, label):
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "true", "t", "yes", "y", "on"}:
            return True
        if normalized in {"0", "false", "f", "no", "n", "off"}:
            return False
    _fail(f"Config key '{label}' must be boolean.")


def _resolve_path(value):
    if value is None or value == "":
        return None
    value = str(value)
    if os.path.isabs(value):
        return os.path.abspath(value)
    return os.path.abspath(os.path.join(WORKDIR, value))


def _is_unset(value):
    return value is None or str(value).strip().lower() in {"", "na", "none", "null"}


def _resolve_optional_path(value):
    if _is_unset(value):
        return None
    return _resolve_path(value)


def _join_path(base, *parts):
    normalized_parts = [str(base)]
    for part in parts:
        normalized_parts.extend(
            piece
            for piece in str(part).replace("\\", "/").split("/")
            if piece
        )
    return os.path.normpath(os.path.join(*normalized_parts))


def outpath(*parts):
    return _join_path(config["output"], *parts)


def tmp_path(base, *parts):
    return _join_path(base, *parts)


def _ensure_file(path, label):
    if not path or not os.path.exists(path):
        _fail(f"{label}: {path} does not exist.")
    if not os.path.isfile(path):
        _fail(f"{label}: {path} is not a file.")
    if not os.access(path, os.R_OK):
        _fail(f"{label}: {path} is not readable.")


def _ensure_parent_writable(path, label):
    parent = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(parent):
        _fail(f"{label}: parent directory {parent} does not exist.")
    if not os.access(parent, os.W_OK):
        _fail(f"{label}: parent directory {parent} is not writable.")


def _is_fasta(path):
    return path.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"))


def _is_bgzip(path):
    return path.endswith(".gz")


def _validate_fasta(path, label, require_pansn=False, require_unique=False):
    _ensure_file(path, label)
    if not _is_fasta(path):
        _fail(f"{label}: {path} must end with .fa, .fasta, .fna, or the .gz versions.")
    _ensure_file(path + ".fai", f"{label} FASTA index")
    if _is_bgzip(path):
        _ensure_file(path + ".gzi", f"{label} bgzip index")
    names = _read_fai_names(path)
    if require_pansn:
        seen = set()
        for name in names:
            if len(name.split("#")) != 3:
                _fail(f"{label}: contig '{name}' does not follow PanSN sample#hap#contig naming.")
            if require_unique and name in seen:
                _fail(f"{label}: duplicate FASTA contig '{name}'.")
            seen.add(name)
    return names


def _read_fai_names(path):
    names = []
    with open(path + ".fai", "r") as handle:
        for line in handle:
            if not line.strip():
                continue
            names.append(line.rstrip().split("\t")[0])
    return names


def _validate_alignment(path, label):
    _ensure_file(path, label)
    if path.endswith(".bam"):
        if not (os.path.exists(path + ".bai") or os.path.exists(path + ".csi")):
            _fail(f"{label}: {path} is missing .bai or .csi index.")
    elif path.endswith(".cram"):
        if not os.path.exists(path + ".crai"):
            _fail(f"{label}: {path} is missing .crai index.")
    else:
        _fail(f"{label}: {path} must be BAM or CRAM.")


def _validate_flagger(path):
    if path is None:
        return
    _ensure_file(path, "Flagger blacklist")
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            contig = line.split("\t", 1)[0]
            if len(contig.split("#")) != 3:
                _fail(
                    f"Flagger blacklist line {lineno}: contig '{contig}' does not follow PanSN naming."
                )


def _read_tsv(path, required_columns, schema_path, label):
    _ensure_file(path, label)
    rows = []
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            _fail(f"{label}: empty TSV.")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            _fail(f"{label}: missing required column(s): {', '.join(missing)}.")
        for lineno, row in enumerate(reader, start=2):
            if row is None or all((value is None or value == "") for value in row.values()):
                continue
            cleaned = {key: (value.strip() if value is not None else value) for key, value in row.items()}
            validate(cleaned, schema_path)
            rows.append(cleaned)
    if not rows:
        _fail(f"{label}: no data rows found.")
    return rows


def _parse_regions(path, reference_contigs):
    _ensure_file(path, "Regions BED")
    rows = OrderedDict()
    with open(path, "r") as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                _fail(f"Regions BED line {lineno}: expected at least 3 tab-separated columns.")
            chrom, start, end = fields[0], fields[1], fields[2]
            if chrom not in reference_contigs:
                _fail(f"Regions BED line {lineno}: chromosome '{chrom}' is absent from the reference .fai.")
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                _fail(f"Regions BED line {lineno}: start and end must be integers.")
            if start_i < 0 or end_i <= start_i:
                _fail(f"Regions BED line {lineno}: expected 0 <= start < end.")
            annot = fields[3] if len(fields) >= 4 and fields[3] else "unknown"
            alts = None
            ploidy_text = None
            field4 = fields[4] if len(fields) >= 5 else ""
            field5 = fields[5] if len(fields) >= 6 else ""
            if field4 and field4 != ".":
                if len(fields) == 5 and field4.isdigit():
                    ploidy_text = field4
                else:
                    alts = field4
            if field5 and field5 != ".":
                ploidy_text = field5
            if len(fields) > 6:
                _fail(f"Regions BED line {lineno}: expected at most 6 columns.")
            if ploidy_text is None:
                ploidy_text = "2"
            try:
                ploidy_i = int(ploidy_text)
            except ValueError:
                _fail(f"Regions BED line {lineno}: ploidy must be a positive integer.")
            if ploidy_i < 1:
                _fail(f"Regions BED line {lineno}: ploidy must be >= 1.")
            if alts is not None:
                _parse_alt_regions(alts, lineno)
            region = f"{chrom}_{start}_{end}"
            if region in rows:
                _fail(f"Regions BED line {lineno}: duplicate region '{region}'.")
            rows[region] = {
                "chrom": chrom,
                "start": str(start_i),
                "end": str(end_i),
                "annot": annot,
                "alts": alts,
                "ploidy": str(ploidy_i),
            }
    if not rows:
        _fail("Regions BED: no regions found.")
    return rows


def _parse_alt_regions(text, lineno):
    parsed = []
    for alt in text.split(","):
        alt = alt.strip()
        if not alt:
            continue
        last_colon = alt.rfind(":")
        if last_colon < 1:
            _fail(f"Regions BED line {lineno}: invalid alt region '{alt}', expected chrom:start-end.")
        alt_chrom = alt[:last_colon]
        coords = alt[last_colon + 1 :]
        if "-" not in coords:
            _fail(f"Regions BED line {lineno}: invalid alt region '{alt}', expected chrom:start-end.")
        start, end = coords.split("-", 1)
        try:
            start_i = int(start)
            end_i = int(end)
        except ValueError:
            _fail(f"Regions BED line {lineno}: invalid alt coordinates in '{alt}'.")
        if start_i < 0 or end_i <= start_i:
            _fail(f"Regions BED line {lineno}: invalid alt coordinates in '{alt}'.")
        parsed.append((alt_chrom, str(start_i), str(end_i)))
    return parsed


def _find_optimal_bindings(paths, min_coverage_threshold=2):
    paths = [path for path in paths if path]
    if not paths:
        return []
    dir_coverage = {}
    for path in paths:
        current_path = ""
        for component in path.split("/"):
            if not component:
                continue
            current_path = f"{current_path}/{component}" if current_path else f"/{component}"
            dir_coverage[current_path] = dir_coverage.get(current_path, 0) + 1
    bindings = set()
    for path in paths:
        best_binding = path
        current_path = ""
        for component in path.split("/"):
            if not component:
                continue
            current_path = f"{current_path}/{component}" if current_path else f"/{component}"
            if dir_coverage[current_path] >= min_coverage_threshold:
                best_binding = current_path
        bindings.add(best_binding)
    optimized = set()
    for binding in sorted(bindings, key=len):
        if any(binding.startswith(existing + "/") for existing in optimized):
            continue
        optimized.add(binding)
    if "/" in optimized and len(optimized) > 1:
        optimized.remove("/")
    return sorted(optimized, key=len)


def _append_container_bind_env(bindings):
    for variable in ("APPTAINER_BINDPATH", "SINGULARITY_BINDPATH", "SINGULARITY_BIND"):
        existing = [x for x in os.environ.get(variable, "").split(",") if x]
        merged = []
        for item in existing + bindings:
            if item not in merged:
                merged.append(item)
        if merged:
            os.environ[variable] = ",".join(merged)


def _metadata_region_bed(region):
    row = REGION_ROWS[region]
    return outpath("metadata", "regions", row["chrom"], f"{region}.bed")


def region_chrom(region):
    return REGION_ROWS[region]["chrom"]


def region_ploidy(region):
    return int(REGION_ROWS[region].get("ploidy", "2"))


def truth_graph_path(wildcards):
    """Pre-built graph holding the true haplotypes for this region."""
    try:
        return TRUTH_GRAPHS[wildcards.region]["graph"]
    except KeyError:
        _fail(f"No truth graph configured for region '{wildcards.region}'.")


def region_annotation(region):
    """Label from column 4 of the regions BED, used as gene_name in reports."""
    return REGION_ROWS[region].get("annot", "unknown")


def region_supports_haploid_diploid_downstream(region):
    return region_ploidy(region) <= 2


def region_bed_path(wildcards):
    return outpath("metadata", "regions", wildcards.chr, f"{wildcards.region}.bed")


def sample_alignment_path(wildcards):
    try:
        return SAMPLES[wildcards.sample]["alignment"]
    except KeyError:
        _fail(f"Unknown sample '{wildcards.sample}'.")


def assembly_fasta_path(wildcards):
    try:
        return ASSEMBLIES[wildcards.chr]["fasta"]
    except KeyError:
        _fail(f"No assembly FASTA configured for chromosome '{wildcards.chr}'.")


def assembly_fai_path(wildcards):
    return assembly_fasta_path(wildcards) + ".fai"


def custom_allele_fasta_path(wildcards):
    try:
        return ALLELES[wildcards.region]["fasta"]
    except KeyError:
        _fail(f"No custom allele FASTA configured for region '{wildcards.region}'.")


def realigned_alignment_path(wildcards):
    if READ_MODE == "short":
        folder = "bwa-mem2"
    elif READ_MODE == "ancient":
        folder = "bwa"
    elif LONG_READ_PRESET is not None:
        folder = f"minimap2/{READ_MODE_LABEL}"
    else:
        _fail(f"Unsupported read_mode '{READ_MODE}'.")
    return outpath(folder, wildcards.sample, wildcards.chr, wildcards.region, f"{wildcards.region}.realigned.cram")


def long_read_preset():
    if LONG_READ_PRESET is None:
        _fail(f"read_mode '{READ_MODE}' does not select a minimap2 long-read preset.")
    return LONG_READ_PRESET


def _read_mode_label(read_mode):
    return read_mode.replace(":", "_").replace("/", "_")


def cluster_args():
    """
    Render config['cluster'] as the --key value option string cluster.r expects.
    Booleans become the lowercase true/false spellings the script parses.
    """
    parts = []
    for key in CLUSTER_DEFAULTS:
        value = config["cluster"][key]
        if isinstance(value, bool):
            value = "true" if value else "false"
        parts.append("--{} {}".format(key.replace("_", "-"), shlex.quote(str(value))))
    return " ".join(parts)


_deep_defaults(config, RESOURCE_DEFAULTS)
_deep_defaults(config, {"cluster": CLUSTER_DEFAULTS})
_migrate_time_to_runtime(config)

config.setdefault("pansn_prefix", "grch38#1#")
config.setdefault("tmpdir", "/tmp")
for key, default in BOOLEAN_DEFAULTS.items():
    config[key] = _normalize_bool(config.get(key, default), key)

unknown_cluster_keys = sorted(set(config["cluster"]) - set(CLUSTER_DEFAULTS))
if unknown_cluster_keys:
    _fail(f"Config key 'cluster': unknown option(s): {', '.join(unknown_cluster_keys)}.")
for key in CLUSTER_BOOLEAN_KEYS:
    config["cluster"][key] = _normalize_bool(config["cluster"][key], f"cluster.{key}")

validate(config, CONFIG_SCHEMA)

READ_MODE = config["read_mode"]
ALLELE_SOURCE = config["allele_source"]
# leave_zero_out: the genotyped samples' own haplotypes stay in the graph, so an
#   exact hit is possible and QV measures whether cosigt finds it.
# leave_all_out: they are removed from the genotyping graph, so cosigt must
#   reconstruct each sample from other people's haplotypes. A second graph
#   containing everything supplies the truth sequences for comparison.
BENCHMARK_MODE = config.get("benchmark_mode", "leave_zero_out")
LONG_READ_PRESET = READ_MODE.split(":", 1)[1] if READ_MODE.startswith("long:") else None
READ_MODE_LABEL = _read_mode_label(READ_MODE)

config["output"] = _resolve_path(config["output"])
config["reference"] = _resolve_path(config["reference"])
config["samples_table"] = _resolve_path(config["samples"])
config["regions_bed"] = _resolve_path(config["regions"])
config["assemblies_table"] = _resolve_optional_path(config.get("assemblies"))
config["alleles_table"] = _resolve_optional_path(config.get("alleles"))
config["truth_graphs_table"] = _resolve_optional_path(config.get("truth_graphs"))
config["gtf"] = _resolve_optional_path(config.get("gtf")) or "NA"
config["proteins"] = _resolve_optional_path(config.get("proteins")) or "NA"
config["flagger_source"] = _resolve_optional_path(config.get("flagger_blacklist"))
config["tmpdir"] = _resolve_path(config.get("tmpdir"))
config["pggb"]["tmpdir"] = _resolve_path(config.get("pggb", {}).get("tmpdir") or config["tmpdir"])
config["flagger_blacklist"] = outpath("metadata", "flagger", "flagger_blacklist.bed")
config["all_regions"] = outpath("metadata", "regions", "all_regions.tsv")

_ensure_parent_writable(config["output"], "Output directory")
REFERENCE_CONTIGS = set(_validate_fasta(config["reference"], "Reference FASTA"))

if (config.get("gtf") == "NA") != (config.get("proteins") == "NA"):
    _fail("Config keys 'gtf' and 'proteins' must be provided together or omitted together.")
if config.get("gtf") != "NA":
    _ensure_file(config.get("gtf"), "GTF")
    if not config.get("gtf").endswith((".gtf", ".gtf.gz")):
        _fail("GTF must end with .gtf or .gtf.gz.")
    _validate_fasta(config.get("proteins"), "Proteins FASTA")

_validate_flagger(config.get("flagger_source"))

REGION_ROWS = _parse_regions(config["regions_bed"], REFERENCE_CONTIGS)
REGION_ORDER = list(REGION_ROWS.keys())
CHROMOSOMES = list(OrderedDict((row["chrom"], None) for row in REGION_ROWS.values()).keys())

sample_rows = _read_tsv(
    config["samples_table"],
    ["sample", "alignment"],
    SAMPLES_SCHEMA,
    "Samples TSV",
)
SAMPLES = OrderedDict()
for row in sample_rows:
    sample = row["sample"]
    if sample in SAMPLES:
        _fail(f"Samples TSV: duplicate sample '{sample}'.")
    alignment = _resolve_path(row["alignment"])
    _validate_alignment(alignment, f"Sample '{sample}' alignment")
    SAMPLES[sample] = {"alignment": alignment}

ASSEMBLIES = OrderedDict()
ALLELES = OrderedDict()
if ALLELE_SOURCE == "assemblies":
    if config.get("assemblies_table") is None:
        _fail("Config key 'assemblies' is required when allele_source is 'assemblies'.")
    assembly_rows = _read_tsv(
        config["assemblies_table"],
        ["chromosome", "fasta"],
        ASSEMBLIES_SCHEMA,
        "Assemblies TSV",
    )
    for row in assembly_rows:
        chrom = row["chromosome"]
        if chrom in ASSEMBLIES:
            _fail(f"Assemblies TSV: duplicate chromosome '{chrom}'.")
        fasta = _resolve_path(row["fasta"])
        _validate_fasta(fasta, f"Assembly FASTA for '{chrom}'", require_pansn=True)
        ASSEMBLIES[chrom] = {"fasta": fasta}
    missing = [chrom for chrom in CHROMOSOMES if chrom not in ASSEMBLIES]
    if missing:
        _fail(f"Assemblies TSV: missing chromosome(s) used by regions: {', '.join(missing)}.")
else:
    if config.get("alleles_table") is None:
        _fail("Config key 'alleles' is required when allele_source is 'custom'.")
    allele_rows = _read_tsv(
        config["alleles_table"],
        ["region", "fasta"],
        ALLELES_SCHEMA,
        "Alleles TSV",
    )
    for row in allele_rows:
        region = row["region"]
        if region in ALLELES:
            _fail(f"Alleles TSV: duplicate region '{region}'.")
        fasta = _resolve_path(row["fasta"])
        _validate_fasta(fasta, f"Custom allele FASTA for '{region}'", require_pansn=True, require_unique=True)
        ALLELES[region] = {"fasta": fasta}
    missing = [region for region in REGION_ORDER if region not in ALLELES]
    if missing:
        _fail(f"Alleles TSV: missing region(s) used by regions BED: {', '.join(missing)}.")

# leave_all_out compares against haplotypes that are, by design, absent from the
# genotyping graph. Those come from pre-built per-region truth graphs supplied by
# the user, one .og per region.
TRUTH_GRAPHS = OrderedDict()
if BENCHMARK_MODE == "leave_all_out":
    if config.get("truth_graphs_table") is None:
        _fail(
            "Config key 'truth_graphs' is required when benchmark_mode is "
            "'leave_all_out': it maps each region to a graph containing the "
            "genotyped samples' own haplotypes, which the genotyping graph "
            "deliberately lacks."
        )
    truth_rows = _read_tsv(
        config["truth_graphs_table"],
        ["region", "graph"],
        TRUTH_GRAPHS_SCHEMA,
        "Truth graphs TSV",
    )
    for row in truth_rows:
        region = row["region"]
        if region in TRUTH_GRAPHS:
            _fail(f"Truth graphs TSV: duplicate region '{region}'.")
        graph = _resolve_path(row["graph"])
        _ensure_file(graph, f"Truth graph for '{region}'")
        if not graph.endswith(".og"):
            _fail(
                f"Truth graph for '{region}': {graph} must be an odgi graph (.og). "
                "Convert a GFA with 'odgi build -g in.gfa -o out.og'."
            )
        TRUTH_GRAPHS[region] = {"graph": graph}
    missing = [region for region in REGION_ORDER if region not in TRUTH_GRAPHS]
    if missing:
        _fail(
            "Truth graphs TSV: missing region(s) used by the regions BED: "
            + ", ".join(missing)
            + "."
        )

config["samples"] = list(SAMPLES.keys())
config["regions"] = REGION_ORDER
config["chromosomes"] = CHROMOSOMES

REGION_BED_TARGETS = [_metadata_region_bed(region) for region in REGION_ORDER]
METADATA_TARGETS = [config["all_regions"], config["flagger_blacklist"]] + REGION_BED_TARGETS
# APPTAINER_ARGS_FILE is appended below, once the helpers it needs are defined.

bind_paths = {
    os.path.dirname(config["reference"]),
    os.path.dirname(config["output"]),
    config["tmpdir"],
    config["pggb"]["tmpdir"],
}
if config.get("gtf") != "NA":
    bind_paths.add(os.path.dirname(config.get("gtf")))
    bind_paths.add(os.path.dirname(config.get("proteins")))
if config.get("flagger_source") is not None:
    bind_paths.add(os.path.dirname(config.get("flagger_source")))
for sample in SAMPLES.values():
    bind_paths.add(os.path.dirname(sample["alignment"]))
for assembly in ASSEMBLIES.values():
    bind_paths.add(os.path.dirname(assembly["fasta"]))
for allele in ALLELES.values():
    bind_paths.add(os.path.dirname(allele["fasta"]))
for truth in TRUTH_GRAPHS.values():
    bind_paths.add(os.path.dirname(truth["graph"]))

BIND_PATHS = _find_optimal_bindings(sorted(bind_paths))
_append_container_bind_env(BIND_PATHS)


def deployment_methods():
    """
    Names of the software deployment methods Snakemake was invoked with, e.g.
    {"apptainer"} or {"conda"}. Empty when running with tools taken from PATH.
    """
    try:
        methods = workflow.deployment_settings.deployment_method
    except AttributeError:
        return set()
    return {getattr(method, "name", str(method)).lower() for method in methods}


def apptainer_args():
    """
    Compose the --apptainer-args string: bind mounts covering every configured
    input/output location, plus -e (--cleanenv), which pggb requires. Extra
    flags can be appended through the 'apptainer_extra' config key.
    """
    parts = []
    if BIND_PATHS:
        parts.append("-B " + ",".join(BIND_PATHS))
    if config.get("apptainer_cleanenv", True):
        parts.append("-e")
    extra = str(config.get("apptainer_extra", "") or "").strip()
    if extra:
        parts.append(extra)
    return " ".join(parts)


def required_tools():
    """
    Binaries the current configuration will actually invoke. Kept in the
    workflow rather than the Makefile because it depends on read_mode,
    allele_source and the optional output switches.
    """
    tools = {
        "samtools", "bgzip", "bedtools", "minimap2", "odgi", "impg",
        "gafpack", "gfainject", "panplexity", "pggb", "cosigt", "Rscript",
    }
    if READ_MODE == "short":
        tools.update({"bwa-mem2", "kfilt", "meryl"})
    elif READ_MODE == "ancient":
        tools.add("bwa")
    if config.get("vcf"):
        tools.update({"bcftools", "tabix", "python"})
    if config.get("gtf") != "NA" and config.get("pangene_viz"):
        tools.update({"pangene", "pangene.js", "miniprot"})
    if config.get("wally_viz"):
        tools.add("wally")
    if config.get("sv_calling"):
        tools.add("svim-asm")
    return sorted(tools)


def check_tools_on_path():
    """
    With no container or conda deployment every tool must already be on PATH.
    Checked here, where the config is known, so the Makefile does not have to
    replicate the read_mode/optional-output logic.
    """
    missing = [tool for tool in required_tools() if which(tool) is None]
    if missing:
        _fail(
            "SOFTWARE=none was requested but these tools are not on PATH: "
            + ", ".join(missing)
            + ". Install them, or use SOFTWARE=apptainer or SOFTWARE=conda."
        )


DEPLOYMENT = deployment_methods()
APPTAINER_ARGS_FILE = os.path.join(WORKDIR, ".cosigt", "apptainer.args")
METADATA_TARGETS.append(APPTAINER_ARGS_FILE)

# Fail fast, alongside the other input validation, rather than part-way through
# a run. Only meaningful without containers or conda, where tools come from PATH.
if not DEPLOYMENT:
    check_tools_on_path()
