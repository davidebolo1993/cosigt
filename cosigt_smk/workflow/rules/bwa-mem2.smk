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
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa-mem2 index {input}
		'''

rule bwa_mem2_samtools_sort:
	'''
	bwa mem and sam-to-bam conversion with samtools
	'''
	input:
		ref=rules.odgi_paths.output,
		idx=rules.bwa_mem2_index.output,
		sample=rules.samtools_fasta.output
	output:
		'results/cosigt_results/{sample}/{sample}.region.realigned.bam'
	threads:
		10
	resources:
		mem_mb=10000
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa-mem2 mem -t {threads} {input.ref} {input.sample} | samtools sort > {output}
		'''