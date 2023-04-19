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

#rule bwa_mem_samtools_sort:
#	'''
#	Bwa mem and bam conversion
#	'''
#	input:
#		ref=rules.odgi_paths.output,
#		idx=rules.bwa_index.output,
#		sample=rules.samtools_fasta.output
#	output:
#		"results/cosigt_results/{sample}/{sample}.region.realigned.bam"
#	threads:
#		10
#	container:
#		'docker://davidebolo1993/graph_genotyper:latest'
#	shell:
#		'''
#		bwa mem -t {threads} {input.ref} {input.sample} | samtools sort > {output}
#		'''

#rule bwa_mem_samtools_sort_old:
#	'''
#	Bwa mem and bam conversion
#	'''
#	input:
#		ref=rules.odgi_paths.output,
#		idx=rules.bwa_index.output,
#		sample=rules.samtools_fasta.output
#	output:
#		"results/cosigt_results/{sample}/{sample}.region.realigned.bam"
#	threads:
#		10
#	container:
#		'docker://davidebolo1993/graph_genotyper:latest'
#	shell:
#		'''
#		bwa mem -t {threads} {input.ref} {input.sample} -k 19 -r 2.5 | samtools sort > {output}
#		'''

rule bwa_aln_old:
	'''
	Bwa aln
	'''
	input:
		ref=rules.odgi_paths.output,
		idx=rules.bwa_index.output,
		sample=rules.samtools_fasta.output
	output:
		"results/cosigt_results/{sample}/{sample}.region.realigned.sai"
	threads:
		10
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa aln -0 -t {threads} -l 1024 -n 0.01 -o 2 {input.ref} {input.sample} > {output}
		'''

rule bwa_samse_old:
	'''
	Bwa samse and samtools sort
	'''
	input:
		sai=rules.bwa_aln_old.output,
		ref=rules.odgi_paths.output,
		sample=rules.samtools_fasta.output
	output:
		"results/cosigt_results/{sample}/{sample}.region.realigned.bam"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa samse {input.ref} {input.sai} {input.sample} | samtools sort > {output}
		'''