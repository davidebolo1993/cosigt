rule bedtools_merge:
	'''
	https://github.com/arq5x/bedtools2
	'''
	input:
		rules.concatenate_batches_per_region.output
	output:
		config['output'] + '/bedtools/{chr}/{region}.bed'
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
		'benchmarks/{chr}.{region}.bedtools_merge.benchmark.txt'
	shell:
		'''
		bedtools sort \
		-i {input} | \
		bedtools merge \
		-d 200000 \
		-i - > {output}
		'''

rule filter_outliers:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.bedtools_merge.output
	output:
		config['output'] + '/bedtools/{chr}/filtered/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.filter_outliers.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/outliers.r \
		{input} \
		{output}
		'''

rule bedtools_getfasta:
	'''
	https://github.com/arq5x/bedtools2
	'''
	input:
		asm_fai=lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
		ref_fasta=rules.pansnspec_target.output,
		ref_fai=rules.samtools_faidx_target.output,
		asm_bed=rules.filter_outliers.output,
		ref_bed=rules.make_reference_bed.output
	output:
		config['output'] + '/samtools/faidx/{chr}/{region}.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	benchmark:
		'benchmarks/{chr}.{region}.bedtools_getfasta.benchmark.txt'
	conda:
		'../envs/bedtools.yaml'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		asm_fai=$(echo {input.asm_fai})
		asm_fasta=$(echo "${{asm_fai%.*}}")
		bedtools getfasta \
			-fi $asm_fasta \
			-bed {input.asm_bed} > {output}
		bedtools getfasta \
			-fi {input.ref_fasta} \
			-bed <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {input.ref_bed}) >> {output}
		'''
