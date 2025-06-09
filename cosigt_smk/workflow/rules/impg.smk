rule impg_project_batches:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.wfmash_align_batches.output,
		bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		unfiltered=config['output'] + '/impg/{chr}/batches/{region}/{batch}.bedpe',
		noblck=config['output'] + '/impg/{chr}/batches/{region}/{batch}.noblck.bedpe',
		merged=config['output'] + '/impg/{chr}/batches/{region}/{batch}.noblck.merged.bedpe',
		filtered=config['output'] + '/impg/{chr}/batches/{region}/{batch}.noblck.merged.filtered.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/impg:0.2.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.{batch}.impg_project_batches.benchmark.txt'
	params:
		blacklist=config['blacklist'],
		flagger_blacklist=config['flagger_blacklist'],
		pansn=config['pansn_prefix'],
		region='{region}'
	shell:
		'''
		impg \
		query \
		-p {input.paf} \
		-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {input.bed}) > {output.unfiltered}
		(grep -v \
		-E \
		-f {params.blacklist} {output.unfiltered} | bedtools \
		intersect \
		-a - \
		-b {params.flagger_blacklist} \
		-F 0.20 \
		-v \
		-wa | bedtools sort -i - || true) > {output.noblck}
		(sh workflow/scripts/bedpe_merge.sh {output.noblck} 200000 || true) > {output.merged}
		(sh workflow/scripts/bedpe_filter.sh {output.merged} {params.region} 2000 || true) > {output.filtered}
		'''

rule concatenate_batches_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/impg/{chr}/batches/{region}/{batch}.noblck.merged.filtered.bedpe', chr='{chr}', region='{region}', batch=get_batches(wildcards))
	output:
		config['output'] + '/impg/{chr}/merged/{region}.bedpe'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.concatenate_batches_per_region.txt'
	shell:
		'''
		cat {input} > {output}
		'''

