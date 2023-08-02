import os

rule odgi_build:
	'''
	odgi build
	'''
	input:
		config['graph'],
	output:
		'results/odgi/build/' + os.path.basename(config['graph']).replace('.gfa', '.og')
	threads:
		config['odgi']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=config['odgi']['mem_mb'],
		time=config['odgi']['time']
	shell:
		'''
		odgi build \
			-g {input} \
			-o {output} \
			-t {threads}
		'''

rule odgi_extract:
	'''
	odgi extract
	'''
	input:
		graph=rules.odgi_build.output,
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		'results/odgi/extract/{region}.og'
	threads:
		config['odgi']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=config['odgi']['mem_mb'],
		time=config['odgi']['time']
	shell:
		'''
		odgi extract \
			-i {input.graph} \
			-o {output} \
			-b {input.region} \
			-d 10000 \
			-t {threads}
		'''
	
rule odgi_chop:
	'''
	odgi chop
	'''
	input:
		rules.odgi_extract.output
	output:
		'results/odgi/chop/{region}.og'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi chop  \
		-i {input} \
		-c 32 \
		-o {output}
		'''

rule odgi_view:
	'''
	odgi view
	'''
	input:
		rules.odgi_chop.output
	output:
		'results/odgi/view/{region}.gfa'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi view \
		-i {input} \
		-g > {output}
		'''

rule odgi_paths_matrix:
	'''
	odgi build non-binary haplotype matrix
	'''
	input:
		rules.odgi_chop.output
	output:
		'results/odgi/paths/matrix/{region}.tsv.gz'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi paths \
		-i {input} \
		-H | cut -f 1,4- | pigz > {output}
		'''

rule odgi_paths:
	'''
	odgi paths
	'''
	input:
		rules.odgi_extract.output
	output:
		'results/odgi/paths/fasta/{region}.fasta'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi paths \
		-i {input} \
		-f > {output}
		'''

rule odgi_similarity:
	'''
	odgi similarity
	'''
	input:
		rules.odgi_extract.output
	output:
		'results/odgi/similarity/{region}.tsv'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		odgi similarity \
		-i {input} > {output}
		'''	