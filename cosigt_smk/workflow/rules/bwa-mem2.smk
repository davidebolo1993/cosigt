rule bwa_mem2_index:
	'''
	bwa index
	'''
	input:
		rules.pgrtk_filter.output
	output:
		multiext('results/pgrtk/fasta/{region}.fa', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')	
	threads:
		1
	resources:
		mem_mb=config['bwa-mem2']['mem_mb'],
		time=config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa-mem2 index {input}
		'''

rule bwa_mem2_samtools_sort:
	'''
	bwa-mem2 and sam-to-bam conversion with samtools
	'''
	input:
		ref=rules.pgrtk_filter.output,
		idx=rules.bwa_mem2_index.output,
		sample=rules.samtools_fasta.output
	output:
		'results/bwa-mem2/{sample}/{region}.realigned.bam'
	threads:
		config['bwa-mem2']['threads']
	resources:
		mem_mb=config['bwa-mem2']['mem_mb'],
		time=config['bwa-mem2']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		samtools_threads=config['samtools']['threads']
	shell:
		'''
		bwa-mem2 mem -t {threads} {input.ref} {input.sample} | samtools sort -@ {params.samtools_threads} - > {output}
		'''
