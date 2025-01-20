rule make_reference_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.impg_project.output
	output:
		config['output'] + '/samtools/bed/{region}.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{region}.make_reference_bed.benchmark.txt'
	params:
		path=config['path']
	shell:
		'''
		grep {params.path} \
		{input} | \
		cut -f 4-6 | \
		bedtools sort -i - | \
		bedtools merge -i - | \
		sed 's/#/\t/' | \
		rev | \
		cut -f 1-3 | \
		rev > {output}
		'''
	
rule samtools_fasta:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_reference_bed.output
	output:
		config['output'] + '/samtools/fasta/{sample}/{region}.fasta.gz'
	threads:
		config['samtools']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'	
	benchmark:
		'benchmarks/{sample}.{region}.samtools_fasta.benchmark.txt'
	params:
		ref=config['reference']
	shell:
		'''
		samtools view \
		-T {params.ref} \
		-@ {threads} \
		-L {input.bed} \
		-M \
		-b \
		{input.sample} | \
		samtools sort -n | \
		samtools fasta \
		-@ {threads} \
		- | gzip > {output}
		'''

rule samtools_faidx_extract:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		fasta=rules.add_target_to_queries.output,
		bed=rules.bedtools_merge.output
	output:
		config['output'] + '/samtools/faidx/{region}.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	benchmark:
		'benchmarks/{region}.samtools_faidx_extract.benchmark.txt'
	conda:
		'../envs/samtools.yaml'
	shell:
		'''
		samtools faidx \
		-r <(awk '{{print $1":"$2+1"-"$3}}' {input.bed}) \
		{input.fasta} > {output}
		'''

rule samtools_faidx_index:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.samtools_faidx_extract.output
	output:
		config['output'] + '/samtools/faidx/{region}.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{region}.samtools_faidx_index.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''