rule minimap2_index:
	'''
	https://github.com/lh3/minimap2
	'''
	input:
		rules.samtools_faidx_extract.output
	output:
		config['output'] + '/samtools/faidx/{region}.fasta.mmi'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['minimap2']['time']
	container:
		'docker://davidebolo1993/minimap2:2.28'
	conda:
		'../envs/minimap2.yaml'
	benchmark:
		'benchmarks/{region}.minimap2_index.benchmark.txt'
	shell:
		'''
		minimap2 -d {output} {input}
		'''

rule minimap2_samtools_sort:
	'''
	https://github.com/lh3/minimap2
	https://github.com/samtools/samtools
	'''
	input:
		ref_index=rules.minimap2_index.output,
		fasta_sample=rules.samtools_fasta.output
	output:
		config['output'] + '/minimap2/{sample}/{region}.realigned.bam'
	threads:
		config['minimap2']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['minimap2']['time']
	container:
		'docker://davidebolo1993/minimap2:2.28'
	conda:
		'../envs/minimap2.yaml'
	benchmark:
		'benchmarks/{sample}.{region}.minimap2_samtools_sort.benchmark.txt'
	params:
		preset=config['minimap2']['preset']
	shell:
		'''
		minimap2 \
		-a \
		-x {params.preset} \
		-t @{threads} \
		{input.ref_index} \
		{input.fasta_sample} | samtools sort \
		-@ {threads} \
		- > {output}
		'''

