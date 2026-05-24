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


rule write_region_bed:
	'''
	Write one per-region BED, including optional alternative intervals.
	'''
	output:
		outpath("metadata/regions/{chr}/{region}.bed")
	run:
		if wildcards.region not in REGION_ROWS:
			raise WorkflowError(f"Unknown region '{wildcards.region}'.")
		row = REGION_ROWS[wildcards.region]
		if row['chrom'] != wildcards.chr:
			raise WorkflowError(
				f"Region '{wildcards.region}' belongs to chromosome '{row['chrom']}', not '{wildcards.chr}'."
			)
		os.makedirs(os.path.dirname(output[0]), exist_ok=True)
		with open(output[0], 'w') as handle:
			handle.write('\t'.join([row['chrom'], row['start'], row['end'], row['annot']]) + '\n')
			if row['alts'] is not None:
				for alt_chrom, alt_start, alt_end in _parse_alt_regions(row['alts'], wildcards.region):
					handle.write('\t'.join([alt_chrom, alt_start, alt_end, alt_chrom]) + '\n')


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
