rule samtools_fasta_mapped:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_reference_bed.output,
		fasta=config['reference']
	output:
		config['output'] + '/samtools/fasta/{sample}/{chr}/{region}.mapped.fasta.gz'
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
		'benchmarks/{sample}.{chr}.{region}.samtools_fasta_mapped.benchmark.txt'
	shell:
		'''
		samtools view \
		-T {input.fasta} \
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

rule samtools_faidx_index:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.bedtools_getfasta.output
	output:
		config['output'] + '/bedtools/getfasta/{chr}/{region}.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.samtools_faidx_index.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''