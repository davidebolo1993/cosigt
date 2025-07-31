rule gfainject_inject:
	'''
	https://github.com/ekg/gfainject
	- Project the reads-to-contigs alignment into the graph
	'''
	input:
		gfa=rules.odgi_view.output,
		cram=rules.bwamem2_mem_samtools_sort.output
	output:
		config['output'] + '/gfainject/{sample}/{chr}/{region}/{region}.gaf.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/gfainject:0.2.1'
	conda:
		'../envs/gfainject.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.gfainject_inject.benchmark.txt'
	shell:
		'''
		samtools view {input.cram} |
		gfainject \
		--gfa {input.gfa} \
		--sam - \
		--alt-hits 10000 | gzip > {output}
		'''
