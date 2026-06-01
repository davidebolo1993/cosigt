CONFIG_FILE ?= .cosigt.mk

COSIGT_DIR ?= cosigt_smk
SNAKEMAKE ?= snakemake
PYTHON ?= python
PROFILE ?= local
SOFTWARE ?= apptainer
TARGET ?= cosigt
CORES ?= 32
SMK_ARGS ?=
APPTAINER_ARGS ?=
CONDA_MIN_VERSION ?= 24.7.1

-include $(CONFIG_FILE)

PROFILE_NAME := $(patsubst profiles/%,%,$(strip $(PROFILE)))
PROFILE_PATH := $(if $(filter none,$(PROFILE_NAME)),,$(if $(filter profiles/%,$(strip $(PROFILE))),$(strip $(PROFILE)),profiles/$(PROFILE_NAME)))
PROFILE_FLAG := $(if $(strip $(PROFILE_PATH)),--profile $(PROFILE_PATH),)
SOFTWARE_NAME := $(strip $(SOFTWARE))
SOFTWARE_FLAG := $(if $(SOFTWARE_NAME),$(if $(filter none,$(SOFTWARE_NAME)),,--software-deployment-method $(SOFTWARE_NAME)),)
APPTAINER_ARG_FLAGS := $(if $(strip $(APPTAINER_ARGS)),--apptainer-args "$(APPTAINER_ARGS)",)

.PHONY: init check check-dryrun dryrun run check-env check-snakemake check-profile check-software \
	check-slurm-plugin check-lsf-plugin check-cluster-generic-plugin

define REQUIRE_PY_MODULE
@command -v $(PYTHON) >/dev/null 2>&1 || { echo "$(PYTHON) was not found. Activate the environment that contains Snakemake, or set PYTHON=/path/to/python."; exit 1; }; \
$(PYTHON) -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('$(1)') else 1)" || { echo "Missing $(2). Install it with: $(PYTHON) -m pip install $(2)"; exit 1; }
endef

define RUN_SNAKEMAKE
cd $(COSIGT_DIR) && $(SNAKEMAKE) $(1) $(PROFILE_FLAG) --cores $(CORES) $(SOFTWARE_FLAG) $(APPTAINER_ARG_FLAGS) $(SMK_ARGS)
endef

init:
	mkdir -p $(COSIGT_DIR)/config
	cp -n $(COSIGT_DIR)/config/config.yaml.example $(COSIGT_DIR)/config/config.yaml
	cp -n $(COSIGT_DIR)/config/samples.tsv.example $(COSIGT_DIR)/config/samples.tsv
	cp -n $(COSIGT_DIR)/config/regions.bed.example $(COSIGT_DIR)/config/regions.bed
	cp -n $(COSIGT_DIR)/config/assemblies.tsv.example $(COSIGT_DIR)/config/assemblies.tsv
	cp -n $(COSIGT_DIR)/config/alleles.tsv.example $(COSIGT_DIR)/config/alleles.tsv
	@printf '%s\n' \
		'# Local COSIGT make defaults. Command-line make variables override these.' \
		'PROFILE := $(PROFILE_NAME)' \
		'SOFTWARE := $(SOFTWARE_NAME)' \
		'TARGET := $(TARGET)' \
		'CORES := $(CORES)' \
		'COSIGT_DIR := $(COSIGT_DIR)' \
		'SNAKEMAKE := $(SNAKEMAKE)' \
		'PYTHON := $(PYTHON)' \
		> $(CONFIG_FILE)
	@echo "Wrote $(CONFIG_FILE)."
	@echo "Next: edit $(COSIGT_DIR)/config/*.yaml/tsv as needed, then run 'make dryrun'."

check-env: check-snakemake check-profile check-software

check-snakemake:
	@command -v $(SNAKEMAKE) >/dev/null 2>&1 || { echo "$(SNAKEMAKE) was not found in PATH. Activate the Snakemake environment or set SNAKEMAKE=/path/to/snakemake."; exit 1; }

check-profile:
	@if [ "$(PROFILE_NAME)" = "slurm" ]; then \
		$(MAKE) --no-print-directory check-slurm-plugin; \
	elif [ "$(PROFILE_NAME)" = "lsf" ]; then \
		$(MAKE) --no-print-directory check-lsf-plugin; \
	elif [ "$(PROFILE_NAME)" = "cluster-generic" ]; then \
		$(MAKE) --no-print-directory check-cluster-generic-plugin; \
	elif [ "$(PROFILE_NAME)" = "local" ] || [ "$(PROFILE_NAME)" = "none" ] || [ -z "$(PROFILE_NAME)" ]; then \
		:; \
	elif [ ! -d "$(COSIGT_DIR)/$(PROFILE_PATH)" ]; then \
		echo "Profile not found: $(COSIGT_DIR)/$(PROFILE_PATH)"; \
		exit 1; \
	fi

check-software:
	@if printf '%s\n' "$(SOFTWARE_NAME)" | tr ',' ' ' | grep -qw apptainer; then \
		command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1 || { echo "SOFTWARE=apptainer was requested, but neither apptainer nor singularity was found in PATH."; exit 1; }; \
	fi
	@if printf '%s\n' "$(SOFTWARE_NAME)" | tr ',' ' ' | grep -qw conda; then \
		command -v conda >/dev/null 2>&1 || { echo "SOFTWARE=conda was requested, but conda was not found in PATH."; exit 1; }; \
		command -v $(PYTHON) >/dev/null 2>&1 || { echo "$(PYTHON) was not found. Activate the environment that contains Snakemake, or set PYTHON=/path/to/python."; exit 1; }; \
		version=$$(conda --version | awk '{print $$2}'); \
		$(PYTHON) -c 'import re, sys; parse=lambda v: tuple(map(int, (re.findall(r"\d+", v) + ["0", "0", "0"])[:3])); sys.exit(0 if parse(sys.argv[1]) >= parse(sys.argv[2]) else 1)' "$$version" "$(CONDA_MIN_VERSION)" || { echo "Snakemake requires conda >= $(CONDA_MIN_VERSION) when SOFTWARE=conda; found $$version."; exit 1; }; \
	fi

check-slurm-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_slurm,snakemake-executor-plugin-slurm)

check-lsf-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_lsf,snakemake-executor-plugin-lsf)

check-cluster-generic-plugin:
	$(call REQUIRE_PY_MODULE,snakemake_executor_plugin_cluster_generic,snakemake-executor-plugin-cluster-generic)

check: check-env
	$(call RUN_SNAKEMAKE,check)

check-dryrun: check-env
	$(call RUN_SNAKEMAKE,check --dry-run)

dryrun: check-env
	$(call RUN_SNAKEMAKE,$(TARGET) --dry-run)

run: check-env
	$(call RUN_SNAKEMAKE,$(TARGET))
