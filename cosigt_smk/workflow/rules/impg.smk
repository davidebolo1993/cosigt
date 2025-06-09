rule impg_project_batches:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.wfmash_align_batches.output,
		bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		unfiltered=config['output'] + '/impg/{chr}/batches/{region}/{batch}.bedpe',
		filtered=config['output'] + '/impg/{chr}/batches/{region}/{batch}.filtered.bedpe'
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
		(impg \
		query \
		-p {input.paf} \
		-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {input.bed}) | \
		grep -v \
		-E \
		-f {params.blacklist} | bedtools \
		intersect \
		-a - \
		-b {params.flagger_blacklist} \
		-F 0.20 \
		-v \
		-wa || true) > {output.unfiltered}
		sh workflow/scripts/check_flanks.sh {output.unfiltered} {params.region} 0.1 > {output.filtered}
		'''

rule concatenate_batches_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/impg/{chr}/batches/{region}/{batch}.filtered.bedpe', chr='{chr}', region='{region}', batch=get_batches(wildcards))
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

rule make_reference_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.concatenate_batches_per_region.output
	output:
		config['output'] + '/samtools/bed/{chr}/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_reference_bed.benchmark.txt'
	params:
		pansn=config['pansn_prefix'] + '{chr}'
	shell:
		'''
		grep {params.pansn} \
		{input} | \
		cut -f 4-6 | \
		bedtools sort -i - | \
		bedtools merge -i - | \
		sed 's/#/\t/g' | \
		rev | \
		cut -f 1-3 | \
		rev > {output}
		'''
