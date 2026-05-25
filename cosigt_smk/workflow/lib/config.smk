import csv
import os
from collections import OrderedDict
from shutil import copyfile

from snakemake.exceptions import WorkflowError
from snakemake.utils import validate


WORKDIR = os.getcwd()
SCHEMA_DIR = os.path.join(WORKDIR, "workflow", "schemas")
CONFIG_SCHEMA = os.path.join(SCHEMA_DIR, "config.schema.yaml")
SAMPLES_SCHEMA = os.path.join(SCHEMA_DIR, "samples.schema.yaml")
ASSEMBLIES_SCHEMA = os.path.join(SCHEMA_DIR, "assemblies.schema.yaml")
ALLELES_SCHEMA = os.path.join(SCHEMA_DIR, "alleles.schema.yaml")


RESOURCE_DEFAULTS = {
    "bwa-mem2": {"threads": 5, "mem_mb": 10000, "runtime": 4},
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
    "meryl": {"threads": 10, "mem_mb": 25000, "runtime": 20},
    "kfilt": {"threads": 8, "mem_mb": 20000, "runtime": 20},
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
            alts = fields[4] if len(fields) >= 5 and fields[4] else None
            if len(fields) > 5:
                _fail(f"Regions BED line {lineno}: expected at most 5 columns.")
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


_deep_defaults(config, RESOURCE_DEFAULTS)
_migrate_time_to_runtime(config)

config.setdefault("pansn_prefix", "grch38#1#")
config.setdefault("tmpdir", "/tmp")
for key, default in BOOLEAN_DEFAULTS.items():
    config[key] = _normalize_bool(config.get(key, default), key)

validate(config, CONFIG_SCHEMA)

READ_MODE = config["read_mode"]
ALLELE_SOURCE = config["allele_source"]
LONG_READ_PRESET = READ_MODE.split(":", 1)[1] if READ_MODE.startswith("long:") else None
READ_MODE_LABEL = _read_mode_label(READ_MODE)

config["output"] = _resolve_path(config["output"])
config["reference"] = _resolve_path(config["reference"])
config["samples_table"] = _resolve_path(config["samples"])
config["regions_bed"] = _resolve_path(config["regions"])
config["assemblies_table"] = _resolve_optional_path(config.get("assemblies"))
config["alleles_table"] = _resolve_optional_path(config.get("alleles"))
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

config["samples"] = list(SAMPLES.keys())
config["regions"] = REGION_ORDER
config["chromosomes"] = CHROMOSOMES

REGION_BED_TARGETS = [_metadata_region_bed(region) for region in REGION_ORDER]
METADATA_TARGETS = [config["all_regions"], config["flagger_blacklist"]] + REGION_BED_TARGETS

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

BIND_PATHS = _find_optimal_bindings(sorted(bind_paths))
_append_container_bind_env(BIND_PATHS)
