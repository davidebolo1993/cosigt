rule bwa_index:
	'''
	Bwa index
	'''
	input:
		rules.odgi_paths.output
	output:
		multiext("resources/odgi/z.fa", ".bwt", ".pac", ".ann", ".amb", ".sa")	
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa index {input}
		'''

rule bwa_mem_samtools_sort:
	'''
	Bwa mem and bam conversion
	'''
	input:
		ref=rules.odgi_paths.output,
		idx=rules.bwa_index.output,
		sample=rules.samtools_fasta.output
	output:
		"results/cosigt_results/{sample}/{sample}.region.realigned.bam"
	threads:
		10
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa mem -t {threads} {input.ref} {input.sample} | samtools sort > {output}
		'''
