rule samtools_fasta_mapped:
	'''
	https://github.com/samtools/samtools
	- Extract reads aligned to the reference
	- Sort them by name (for smart pairing)
	- Convert to fasta
	- Compress
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_alignment_bed.output,
		fasta=config['reference']
	output:
		temp(config['output'] + '/samtools/fasta/{sample}/{chr}/{region}/{region}.mapped.fasta.gz')
	threads:
		config['samtools']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
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
