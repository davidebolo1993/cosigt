COSIGT_DIR ?= cosigt_smk
SNAKEMAKE ?= snakemake
PYTHON ?= python
PROFILE ?= profiles/slurm
SOFTWARE ?= conda
TARGET ?= cosigt
CORES ?= 32
SMK_ARGS ?=
CONDA_MIN_VERSION ?= 24.7.1

.PHONY: init check dryrun run run-slurm run-lsf run-cluster-generic \
	check-profile-plugin check-slurm-plugin check-lsf-plugin \
	check-cluster-generic-plugin check-conda-version install-cluster-plugins

define REQUIRE_PY_MODULE
@command -v $(PYTHON) >/dev/null 2>&1 || { echo "$(PYTHON) was not found. Activate the environment that contains Snakemake, or set PYTHON=/path/to/python."; exit 1; }; \
$(PYTHON) -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('$(1)') else 1)" || { echo "Missing $(2). Install it with: $(PYTHON) -m pip install $(2)"; exit 1; }
endef

init:
	mkdir -p $(COSIGT_DIR)/config
	cp -n $(COSIGT_DIR)/config/config.yaml.example $(COSIGT_DIR)/config/config.yaml
	cp -n $(COSIGT_DIR)/config/samples.tsv.example $(COSIGT_DIR)/config/samples.tsv
	cp -n $(COSIGT_DIR)/config/regions.bed.example $(COSIGT_DIR)/config/regions.bed
	cp -n $(COSIGT_DIR)/config/assemblies.tsv.example $(COSIGT_DIR)/config/assemblies.tsv
	cp -n $(COSIGT_DIR)/config/alleles.tsv.example $(COSIGT_DIR)/config/alleles.tsv

check-profile-plugin:
	@if [ "$(PROFILE)" = "profiles/slurm" ]; then \
		$(MAKE) --no-print-directory check-slurm-plugin; \
	elif [ "$(PROFILE)" = "profiles/lsf" ]; then \
		$(MAKE) --no-print-directory check-lsf-plugin; \
	elif [ "$(PROFILE)" = "profiles/cluster-generic" ]; then \
		$(MAKE) --no-print-directory check-cluster-generic-plugin; \
	fi

check-slurm-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_slurm,snakemake-executor-plugin-slurm)

check-lsf-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_lsf,snakemake-executor-plugin-lsf)

check-cluster-generic-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_cluster_generic,snakemake-executor-plugin-cluster-generic)

check-conda-version:
	@if printf '%s\n' "$(SOFTWARE)" | tr ',' ' ' | grep -qw conda; then \
		command -v conda >/dev/null 2>&1 || { echo "SOFTWARE=conda was requested, but conda was not found in PATH."; exit 1; }; \
		command -v $(PYTHON) >/dev/null 2>&1 || { echo "$(PYTHON) was not found. Activate the environment that contains Snakemake, or set PYTHON=/path/to/python."; exit 1; }; \
		version=$$(conda --version | awk '{print $$2}'); \
		$(PYTHON) -c 'import re, sys; parse=lambda v: tuple(map(int, (re.findall(r"\d+", v) + ["0", "0", "0"])[:3])); sys.exit(0 if parse(sys.argv[1]) >= parse(sys.argv[2]) else 1)' "$$version" "$(CONDA_MIN_VERSION)" || { echo "Snakemake requires conda >= $(CONDA_MIN_VERSION) when SOFTWARE=conda; found $$version."; echo "Update conda in the Snakemake environment, or run with another software deployment method if available."; exit 1; }; \
	fi

install-cluster-plugins:
	$(PYTHON) -m pip install snakemake-executor-plugin-slurm snakemake-executor-plugin-lsf snakemake-executor-plugin-cluster-generic

check: check-profile-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) --dry-run check $(SMK_ARGS)

dryrun: check-profile-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) --dry-run $(TARGET) $(SMK_ARGS)

run: check-profile-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-slurm: check-slurm-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/slurm --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-lsf: check-lsf-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/lsf --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-cluster-generic: check-cluster-generic-plugin check-conda-version
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/cluster-generic --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)
