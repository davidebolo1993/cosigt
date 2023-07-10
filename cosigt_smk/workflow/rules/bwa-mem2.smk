rule bwa_mem2_index:
	'''
	bwa index
	'''
	input:
		rules.odgi_paths.output
	output:
		multiext('resources/odgi/z.fa', '.bwt.2bit.64', '.pac', '.ann', '.amb', '.0123')	
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
		ref=rules.odgi_paths.output,
		idx=rules.bwa_mem2_index.output,
		sample=rules.samtools_fasta.output
	output:
		'results/cosigt_results/{sample}/{sample}.region.realigned.bam'
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