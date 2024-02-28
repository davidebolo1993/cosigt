rule minimap2_index:
	'''
	minimap2 index
	'''
	input:
		rules.odgi_paths_fasta.output
	output:
		'results/odgi/paths/fasta/{region}.fa.mmi'
	threads:
		1
	resources:
		mem_mb=config['minimap2']['mem_mb'],
		time=config['minimap2']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		minimap2 -d {output} {input}
		'''

rule minimap2_samtools_sort:
	'''
	minimap2 and sam-to-bam conversion with samtools
	'''
	input:
		idx=rules.minimap2_index.output,
		sample=rules.samtools_fasta.output
	output:
		'results/minimap2/{sample}/{region}.realigned.bam'
	threads:
		config['minimap2']['threads']
	resources:
		mem_mb=config['minimap2']['mem_mb'],
		time=config['minimap2']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		samtools_threads=config['samtools']['threads'],
		preset=config['minimap2']['preset']
	shell:
		'''
		minimap2 -a -x {params.preset} -t {threads} {input.idx} {input.sample} | samtools sort -@ {params.samtools_threads} - > {output}
		'''
