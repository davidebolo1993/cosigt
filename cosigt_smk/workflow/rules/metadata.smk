rule write_all_regions:
	'''
	Write the normalized regions table consumed by reporting/refinement rules.
	'''
	output:
		config['all_regions']
	run:
		os.makedirs(os.path.dirname(output[0]), exist_ok=True)
		with open(output[0], 'w') as handle:
			for region in REGION_ORDER:
				row = REGION_ROWS[region]
				handle.write('\t'.join([row['chrom'], row['start'], row['end'], row['annot']]) + '\n')


rule write_region_beds:
	'''
	Write every per-region BED, including optional alternative intervals.

	One job for all regions rather than one job per region. These are `run:`
	directives, and Snakemake re-establishes the workflow context for each such
	job -- around a second apiece, against a few milliseconds of actual work --
	so a per-region rule made `make check` scale with the size of the regions
	BED and dominate its runtime. Nothing downstream consumes these files
	individually, so there is nothing to gain from the finer granularity.
	'''
	output:
		REGION_BED_TARGETS
	run:
		for region in REGION_ORDER:
			row = REGION_ROWS[region]
			path = _metadata_region_bed(region)
			os.makedirs(os.path.dirname(path), exist_ok=True)
			with open(path, 'w') as handle:
				handle.write('\t'.join([row['chrom'], row['start'], row['end'], row['annot']]) + '\n')
				if row['alts'] is not None:
					for alt_chrom, alt_start, alt_end in _parse_alt_regions(row['alts'], region):
						handle.write('\t'.join([alt_chrom, alt_start, alt_end, alt_chrom]) + '\n')


rule write_apptainer_args:
	'''
	Compose the Apptainer/Singularity flags for this configuration and write
	them where the Makefile can pick them up. Bind mounts are derived from every
	configured input and output location, so users do not have to work them out
	by hand; -e (--cleanenv) is included because pggb fails without it.
	Tools required from PATH when running without containers or conda are
	verified at config-parse time, in lib/config.smk.

	The config and the tables it points to are declared as inputs so that
	editing them rebuilds the flags. Without that this rule has neither input
	nor params, so nothing can invalidate its output and the flags silently
	keep whatever paths were configured the first time `check` ran.
	'''
	input:
		APPTAINER_ARGS_DEPS
	output:
		APPTAINER_ARGS_FILE
	run:
		os.makedirs(os.path.dirname(output[0]), exist_ok=True)
		with open(output[0], 'w') as handle:
			handle.write(apptainer_args() + '\n')


rule write_flagger_blacklist:
	'''
	Write an empty or copied flagger blacklist at a workflow-owned path.
	'''
	input:
		lambda wildcards: config.get('flagger_source') or []
	output:
		config['flagger_blacklist']
	run:
		os.makedirs(os.path.dirname(output[0]), exist_ok=True)
		if len(input) == 0:
			open(output[0], 'w').close()
		else:
			copyfile(input[0], output[0])
