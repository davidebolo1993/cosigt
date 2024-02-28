rule bwa_index:
	'''
	bwa index
	'''
	input:
		rules.odgi_paths_fasta.output
	output:
		multiext('results/odgi/paths/fasta/{region}.fa', '.bwt', '.pac', '.ann', '.amb', '.sa')	
	threads:
		1
	resources:
		mem_mb=config['bwa']['mem_mb'],
		time=config['bwa']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa index {input}
		'''

rule bwa_aln_old:
	'''
	bwa aln
	'''
	input:
		ref=rules.odgi_paths_fasta.output,
		idx=rules.bwa_index.output,
		sample=rules.samtools_fasta.output
	output:
		'results/bwa/{sample}/{region}.realigned.sai'
	threads:
		config['bwa']['threads']
	resources:
		mem_mb=config['bwa']['mem_mb'],
		time=config['bwa']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa aln -0 -t {threads} -l 1024 -n 0.01 -o 2 {input.ref} {input.sample} > {output}
		'''

rule bwa_samse_old_samtools_sort:
	'''
	bwa samse and samtools sort
	'''
	input:
		sai=rules.bwa_aln_old.output,
		ref=rules.odgi_paths_fasta.output,
		sample=rules.samtools_fasta.output
	output:
		'results/bwa/{sample}/{region}.realigned.bam'
	threads:
		config['samtools']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		bwa samse {input.ref} {input.sai} {input.sample} | samtools sort -@ {threads} > {output}
		'''