rule bwa_mem2_index:
	'''
	bwa index
	'''
	input:
		rules.odgi_paths_fasta.output
	output:
		multiext(config['output'] + '/odgi/paths/fasta/{region}.fa', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')	
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	shell:
		'''
		bwa-mem2 index {input}
		'''

rule bwa_mem2_samtools_sort:
	'''
	bwa-mem2 and sam-to-bam conversion with samtools
	'''
	input:
		ref=rules.odgi_paths_fasta.output,
		idx=rules.bwa_mem2_index.output,
		sample=rules.samtools_fasta.output
	output:
		config['output'] + '/bwa-mem2/{sample}/{region}.realigned.bam'
	threads:
		config['bwa-mem2']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['bwa-mem2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	shell:
		'''
		bwa-mem2 mem -t {threads} {input.ref} {input.sample} | samtools sort -@ {threads} - > {output}
		'''
