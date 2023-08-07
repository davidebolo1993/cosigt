rule pgrtk_get_seq:
	input:
		agc=config['agc'],
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		'results/pgrtk/{region}.fa'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		padding=config['pgrtk']['padding']
	shell:
		'''
		pgrtk.py {input.agc} {input.region} {params.padding} {output}
		'''