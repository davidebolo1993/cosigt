rule pgrtk_get_seq:
	'''
	run pgrtk api
	'''
	input:
		agc=config['agc'],
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		'results/pgrtk/{region}.fa'
	threads:
		config['pgrtk']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		padding=config['pgrtk']['padding']
	resources:
		mem_mb=config['pgrtk']['mem_mb'],
		time=config['pgrtk']['time']
	shell:
		'''
		python workflow/scripts/pypgrtk.py {input.agc} {input.region} {params.padding} {output}
		'''

rule faidx:
	'''
	samtools faidx
	'''
	input:
		rules.pgrtk_get_seq.output
	output:
		'results/pgrtk/{region}.fa.fai'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		samtools faidx {input}
		'''