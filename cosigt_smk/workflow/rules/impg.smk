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
	- Merge bedpe entries in a 200k range
	- Filter out contigs that do not span target flanks (1k) and remove contigs that align to multiple blocks
	'''
	input:
		paf=get_merged_paf,
		bed=region_bed_path,
		flagger=rules.write_flagger_blacklist.output,
		index=rules.impg_index.output
	output:
		unfiltered=temp(outpath("impg/{chr}/{region}/{region}.bedpe.gz")),
		noblck=temp(outpath("impg/{chr}/{region}/{region}.noblck.bedpe.gz")),
		merged=temp(outpath("impg/{chr}/{region}/{region}.noblck.merged.bedpe.gz")),
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
		chr='{chr}',
		tmp_bed=outpath("impg/{chr}/{region}/{region}.tmp.bed")
	shell:
		'''
		grep -w {params.chr} {input.bed} > {params.tmp_bed}
		impg \
			query \
			-p {input.paf} \
			-i {input.index} \
			-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {params.tmp_bed}) | gzip > {output.unfiltered}
		rm {params.tmp_bed}
		(bedtools \
			intersect \
			-a {output.unfiltered} \
			-b {input.flagger} \
			-v \
			-wa | bedtools sort -i - || true) | gzip > {output.noblck}
		(zcat {output.noblck} | sh workflow/scripts/bedpe_merge.sh - 200000 || true) | gzip > {output.merged}
		(zcat {output.merged} | sh workflow/scripts/bedpe_filter.sh - {params.region} 1000 || true) | gzip > {output.filtered}
		'''
