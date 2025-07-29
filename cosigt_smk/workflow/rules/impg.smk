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
		bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		unfiltered=config['output'] + '/impg/{chr}/{region}/{region}.bedpe.gz',
		noblck=config['output'] + '/impg/{chr}/{region}/{region}.noblck.bedpe.gz',
		merged=config['output'] + '/impg/{chr}/{region}/{region}.noblck.merged.bedpe.gz',
		filtered=config['output'] + '/impg/{chr}/{region}/{region}.noblck.merged.filtered.bedpe.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/impg:0.2.4'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.impg_project_batches.benchmark.txt'
	params:
		blacklist=config['blacklist'],
		flagger_blacklist=config['flagger_blacklist'],
		pansn=config['pansn_prefix'],
		region='{region}',
		chr='{chr}',
		tmp_bed=config['output'] + '/impg/{chr}/{region}/{region}.tmp.bed'
	shell:
		'''
		grep -w {params.chr} {input.bed} > {params.tmp_bed}
		impg \
			query \
			-p {input.paf} \
			-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {params.tmp_bed}) | gzip > {output.unfiltered}
		rm {params.tmp_bed}
		(zgrep -v \
			-E \
			-f {params.blacklist} {output.unfiltered} | bedtools \
			intersect \
			-a - \
			-b {params.flagger_blacklist} \
			-F 0.1 \
			-v \
			-wa | bedtools sort -i - || true) | gzip > {output.noblck}
		(zcat {output.noblck} | sh workflow/scripts/bedpe_merge.sh - 200000 || true) | gzip > {output.merged}
		(zcat {output.merged} | sh workflow/scripts/bedpe_filter.sh - {params.region} 1000 || true) | gzip > {output.filtered}
		'''
