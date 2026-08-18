                                                                       # cosigt
# Two things to run:
#   make check   validate the environment and the configuration, and work out
#                the Apptainer flags this run needs
#   make run     run the pipeline
#
# Everything is controlled by the variables below. Set them on the command line
# (make run PROFILE=slurm) or persist them in .cosigt.mk via `make init`.

CONFIG_FILE ?= .cosigt.mk

COSIGT_DIR ?= cosigt_smk
SNAKEMAKE  ?= snakemake
PYTHON     ?= python
PROFILE    ?= local
SOFTWARE   ?= apptainer
TARGET     ?= cosigt
SMK_ARGS   ?=
CONDA_MIN_VERSION ?= 24.7.1

-include $(CONFIG_FILE)

# Default to every core the machine reports. Matters most for PROFILE=local,
# where this is the actual parallelism; on cluster profiles it only bounds
# local rules, and the profile's own `jobs:` governs submissions.
DETECTED_CORES := $(shell getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
CORES ?= $(DETECTED_CORES)

PROFILE_NAME  := $(patsubst profiles/%,%,$(strip $(PROFILE)))
PROFILE_PATH  := $(if $(filter none,$(PROFILE_NAME)),,$(if $(filter profiles/%,$(strip $(PROFILE))),$(strip $(PROFILE)),profiles/$(PROFILE_NAME)))
PROFILE_FLAG  := $(if $(strip $(PROFILE_PATH)),--profile $(PROFILE_PATH),)
SOFTWARE_NAME := $(strip $(SOFTWARE))
SOFTWARE_FLAG := $(if $(SOFTWARE_NAME),$(if $(filter none,$(SOFTWARE_NAME)),,--software-deployment-method $(SOFTWARE_NAME)),)

# Written by `make check` (rule write_apptainer_args): bind mounts for every
# configured path, plus -e for pggb.
# Absolute, because RUN_SNAKEMAKE cd's into $(COSIGT_DIR) before the shell
# expands the command substitution that reads this file.
ARGS_FILE := $(abspath $(COSIGT_DIR)/.cosigt/apptainer.args)
USES_APPTAINER := $(shell printf '%s\n' "$(SOFTWARE_NAME)" | tr ',' ' ' | grep -qw apptainer && echo yes)

.PHONY: init check run

define RUN_SNAKEMAKE
cd $(COSIGT_DIR) && $(SNAKEMAKE) $(1) $(PROFILE_FLAG) --cores $(CORES) $(SOFTWARE_FLAG) $(2) $(SMK_ARGS)
endef

define REQUIRE_PY_MODULE
$(PYTHON) -c "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('$(1)') else 1)" \
	|| { echo "  FAIL  missing $(2). Install it with: $(PYTHON) -m pip install $(2)"; exit 1; }; \
echo "  ok    executor plugin $(2)"
endef

init:
	@mkdir -p $(COSIGT_DIR)/config
	@cp -n $(COSIGT_DIR)/config/config.yaml.example    $(COSIGT_DIR)/config/config.yaml
	@cp -n $(COSIGT_DIR)/config/samples.tsv.example    $(COSIGT_DIR)/config/samples.tsv
	@cp -n $(COSIGT_DIR)/config/regions.bed.example    $(COSIGT_DIR)/config/regions.bed
	@cp -n $(COSIGT_DIR)/config/assemblies.tsv.example $(COSIGT_DIR)/config/assemblies.tsv
	@cp -n $(COSIGT_DIR)/config/alleles.tsv.example    $(COSIGT_DIR)/config/alleles.tsv
	@printf '%s\n' \
		'# Local cosigt defaults. Command-line variables override these.' \
		'PROFILE := $(PROFILE_NAME)' \
		'SOFTWARE := $(SOFTWARE_NAME)' \
		'TARGET := $(TARGET)' \
		'CORES := $(CORES)' \
		'COSIGT_DIR := $(COSIGT_DIR)' \
		'SNAKEMAKE := $(SNAKEMAKE)' \
		'PYTHON := $(PYTHON)' \
		> $(CONFIG_FILE)
	@echo "Wrote $(CONFIG_FILE) and $(COSIGT_DIR)/config/."
	@echo "Next: edit $(COSIGT_DIR)/config/*, then run 'make check'."

check:
	@echo "cosigt check  (profile=$(PROFILE_NAME) software=$(SOFTWARE_NAME) target=$(TARGET) cores=$(CORES))"
	@echo
	@echo "environment:"
	@command -v $(SNAKEMAKE) >/dev/null 2>&1 \
		|| { echo "  FAIL  $(SNAKEMAKE) not found in PATH. Activate the Snakemake environment or set SNAKEMAKE=/path/to/snakemake."; exit 1; }
	@echo "  ok    snakemake ($$($(SNAKEMAKE) --version 2>/dev/null))"
	@command -v $(PYTHON) >/dev/null 2>&1 \
		|| { echo "  FAIL  $(PYTHON) not found. Activate the environment that contains Snakemake, or set PYTHON=/path/to/python."; exit 1; }
	@# --- executor plugin required by the selected profile ---
	@case "$(PROFILE_NAME)" in \
	  slurm)           $(call REQUIRE_PY_MODULE,snakemake_executor_plugin_slurm,snakemake-executor-plugin-slurm) ;; \
	  lsf)             $(call REQUIRE_PY_MODULE,snakemake_executor_plugin_lsf,snakemake-executor-plugin-lsf) ;; \
	  cluster-generic) $(call REQUIRE_PY_MODULE,snakemake_executor_plugin_cluster_generic,snakemake-executor-plugin-cluster-generic) ;; \
	  local|none|"")   echo "  ok    local execution, $(CORES) of $(DETECTED_CORES) detected cores" ;; \
	  *) if [ ! -d "$(COSIGT_DIR)/$(PROFILE_PATH)" ]; then \
	       echo "  FAIL  profile not found: $(COSIGT_DIR)/$(PROFILE_PATH)"; exit 1; \
	     else echo "  ok    profile $(PROFILE_PATH)"; fi ;; \
	esac
	@# --- software deployment ---
	@if printf '%s\n' "$(SOFTWARE_NAME)" | tr ',' ' ' | grep -qw apptainer; then \
		if command -v apptainer >/dev/null 2>&1; then echo "  ok    apptainer ($$(apptainer --version 2>/dev/null))"; \
		elif command -v singularity >/dev/null 2>&1; then echo "  ok    singularity ($$(singularity --version 2>/dev/null))"; \
		else echo "  FAIL  SOFTWARE=apptainer requested, but neither apptainer nor singularity is in PATH."; exit 1; fi; \
	fi
	@if printf '%s\n' "$(SOFTWARE_NAME)" | tr ',' ' ' | grep -qw conda; then \
		command -v conda >/dev/null 2>&1 || { echo "  FAIL  SOFTWARE=conda requested, but conda is not in PATH."; exit 1; }; \
		version=$$(conda --version | awk '{print $$2}'); \
		$(PYTHON) -c 'import re, sys; parse=lambda v: tuple(map(int, (re.findall(r"\d+", v) + ["0","0","0"])[:3])); sys.exit(0 if parse(sys.argv[1]) >= parse(sys.argv[2]) else 1)' "$$version" "$(CONDA_MIN_VERSION)" \
			|| { echo "  FAIL  Snakemake needs conda >= $(CONDA_MIN_VERSION) for SOFTWARE=conda; found $$version."; exit 1; }; \
		echo "  ok    conda $$version"; \
	fi
	@if [ "$(SOFTWARE_NAME)" = "none" ]; then echo "  ..    no container/conda deployment; required tools are checked against PATH below"; fi
	@echo
	@echo "configuration and inputs:"
	@log=$$(mktemp); \
	( $(call RUN_SNAKEMAKE,check,) ) >/dev/null 2>$$log \
		|| { echo "  FAIL  see below"; echo; sed 's/^/  /' $$log; rm -f $$log; exit 1; }; \
	rm -f $$log
	@echo "  ok    config, sample table, regions, indexes and input files"
	@echo "  ok    region metadata and flagger blacklist written"
ifeq ($(USES_APPTAINER),yes)
	@echo
	@echo "apptainer flags (used automatically by 'make run'):"
	@sed 's/^/  /' $(ARGS_FILE)
endif
	@echo
	@echo "All checks passed. Run the pipeline with: make run"

run:
	@command -v $(SNAKEMAKE) >/dev/null 2>&1 \
		|| { echo "$(SNAKEMAKE) not found in PATH. Activate the Snakemake environment or set SNAKEMAKE=/path/to/snakemake."; exit 1; }
ifeq ($(USES_APPTAINER),yes)
	@test -s $(ARGS_FILE) \
		|| { echo "$(ARGS_FILE) is missing. Run 'make check' first so the Apptainer flags can be composed."; exit 1; }
	$(call RUN_SNAKEMAKE,$(TARGET),--apptainer-args "$$(cat $(ARGS_FILE))")
else
	$(call RUN_SNAKEMAKE,$(TARGET),)
endif
