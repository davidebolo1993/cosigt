COSIGT_DIR ?= cosigt_smk
SNAKEMAKE ?= snakemake
PYTHON ?= python
PROFILE ?= profiles/local
SOFTWARE ?= conda
TARGET ?= cosigt
CORES ?= 32
SMK_ARGS ?=

.PHONY: init check dryrun run run-slurm run-lsf run-cluster-generic \
	check-profile-plugin check-slurm-plugin check-lsf-plugin \
	check-cluster-generic-plugin install-cluster-plugins

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

install-cluster-plugins:
	$(PYTHON) -m pip install snakemake-executor-plugin-slurm snakemake-executor-plugin-lsf snakemake-executor-plugin-cluster-generic

check: check-profile-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) --dry-run check $(SMK_ARGS)

dryrun: check-profile-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) --dry-run $(TARGET) $(SMK_ARGS)

run: check-profile-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile $(PROFILE) --cores $(CORES) --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-slurm: check-slurm-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/slurm --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-lsf: check-lsf-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/lsf --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)

run-cluster-generic: check-cluster-generic-plugin
	cd $(COSIGT_DIR) && $(SNAKEMAKE) --profile profiles/cluster-generic --software-deployment-method $(SOFTWARE) $(TARGET) $(SMK_ARGS)
