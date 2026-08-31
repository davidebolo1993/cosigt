rule impg_index:
	'''
	https://github.com/pangenome/impg
	- Impg refine	
	'''
	input:
		paf=get_merged_paf
	output:
		outpath("impg/{chr}/{chr}.paf.gz.impg")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://davidebolo1993/impg:0.3.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.impg_index.benchmark.txt'
	shell:
		'''
		impg index \
			-p {input} \
			-i {output} \
			-f
		'''

rule impg_project_batches:
	'''
	https://github.com/pangenome/impg
	- Lift-over regions of interest from the target to the queries
	- Remove user-blacklisted contigs and contigs spanning at least 10% of a flagger bad region
	- Merge bedpe entries in a 200k range and keep only contigs spanning the target flanks (1k)
	'''
	input:
		paf=get_merged_paf,
		bed=region_bed_path,
		flagger=rules.write_flagger_blacklist.output,
		index=rules.impg_index.output
	output:
		unfiltered=temp(outpath("impg/{chr}/{region}/{region}.bedpe.gz")),
		noblck=temp(outpath("impg/{chr}/{region}/{region}.noblck.bedpe.gz")),
		region_bed=temp(outpath("impg/{chr}/{region}/{region}.tmp.bed")),
		filtered=temp(outpath("impg/{chr}/{region}/{region}.noblck.merged.filtered.bedpe.gz"))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt *  config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/impg:0.3.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.impg_project_batches.benchmark.txt'
	params:
		pansn=config['pansn_prefix'],
		region='{region}',
		chr='{chr}'
	shell:
		'''
		grep -w {params.chr} {input.bed} > {output.region_bed}
		impg \
			query \
			-p {input.paf} \
			-i {input.index} \
			-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {output.region_bed}) | gzip > {output.unfiltered}
		# No `|| true` here: every legitimately empty case already exits 0, so it
		# only ever masked real failures, leaving a valid-but-empty gzip and a
		# rule that reported success with no alleles for the region.
		bedtools \
			intersect \
			-a {output.unfiltered} \
			-b {input.flagger} \
			-v \
			-wa | bedtools sort -i - | gzip > {output.noblck}
		zcat {output.noblck} | bash workflow/scripts/bedpe_merge_filter.sh - 200000 {params.region} 1000 | gzip > {output.filtered}
		'''
