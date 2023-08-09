#rule pgrtk_get_seq:
#	'''
#	run pgrtk api - deprecating in favour of cli, but keeping this here just in case want to re-introduce at some point
#	'''
#	input:
#		agc=config['agc'],
#		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
#	output:
#		'results/pgrtk/{region}.fa'
#	threads:
#		config['pgrtk']['threads']
#	container:
#		'docker://davidebolo1993/graph_genotyper:latest'
#	params:
#		padding=config['pgrtk']['padding']
#	resources:
#		mem_mb=config['pgrtk']['mem_mb'],
#		time=config['pgrtk']['time']
#	shell:
#		'''
#		python workflow/scripts/pypgrtk.py {input.agc} {input.region} {params.padding} {output}
#		'''

rule pgrtk_fetch_seqs:
	'''
	run pgrtk fetch
	'''
	input:
		agc=config['agc'],
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		'results/pgrtk/fetch/{region}.fa'
	threads:
		32
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=config['pgrtk']['mem_mb'],
		time=config['pgrtk']['time']
	params:
		agc_prefix=config['agc'].replace('.agc','')
	shell:
		'''
		pgr-fetch-seqs {params.agc_prefix} -r {input.region} > {output}
		'''

rule pgrtk_query:
	'''
	run pgrtk query
	'''
	input:
		agc=config['agc'],
		fasta=rules.pgrtk_fetch_seqs.output
	output:
		hit='results/pgrtk/query/{region}.000.hit',
		fasta='results/pgrtk/query/{region}.000.fa'
	threads:
		32
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=config['pgrtk']['mem_mb'],
		time=config['pgrtk']['time']
	params:
		prefix='results/pgrtk/query/{region}',
		agc_prefix=config['agc'].replace('.agc','')
	shell:
		'''
		pgr-query {params.agc_prefix} {input.fasta} {params.prefix}
		'''

rule pgrtk_filter:
	'''
	filter hits
	'''
	input:
		hit=rules.pgrtk_query.output.hit,
		fasta=rules.pgrtk_query.output.fasta
	output:
		'results/pgrtk/fasta/{region}.fa'
	threads:
		1
	conda:
		'../envs/python.yml'
	shell:
		'''
		python workflow/scripts/filter_pgrtk.py {input.fasta} {input.hit} {output}
		'''

rule faidx:
	'''
	samtools faidx
	'''
	input:
		rules.pgrtk_filter.output
	output:
		'results/pgrtk/fasta/{region}.fa.fai'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		samtools faidx {input}
		'''