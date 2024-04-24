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
	benchmark:
		'benchmarks/{region}.make_reference_bed.benchmark.txt'
	params:
		path=config['path']
	shell:
		'''
		grep {params.path} \
		{input} | \
		sed 's/#/\t/' | \
		rev | \
		cut -f 1-3 | \
		rev > {output}
		'''
	
rule samtools_view:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_reference_bed.output
	output:
		config['output'] + '/samtools/view/{sample}/{region}.bam'
	threads:
		config['samtools']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{sample}.{region}.samtools_view.benchmark.txt'
	params:
		ref=config['reference']
	shell:
		'''
		samtools view \
		-O bam \
		-o {output} \
		-T {params.ref} \
		-@ {threads} \
		-L {input.bed} \
		-M \
		{input.sample}
		'''

rule samtools_fasta:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.samtools_view.output
	output:
		config['output'] + '/samtools/fasta/{sample}/{region}.fasta.gz'
	threads:
		config['samtools']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	benchmark:
		'benchmarks/{sample}.{region}.samtools_fasta.benchmark.txt'
	shell:
		'''
		samtools fasta \
		-@ {threads} \
		{input} | pigz > {output}
		'''