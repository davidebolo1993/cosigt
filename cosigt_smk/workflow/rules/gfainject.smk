rule gfainject_inject:
	'''
	https://github.com/ekg/gfainject
	- Project the reads-to-contigs alignment into the graph
	'''
	input:
		gfa=rules.odgi_utils.output.gfa,
		cram=realigned_alignment_path
	output:
		temp(outpath("gfainject/{sample}/{chr}/{region}/{region}.gaf.gz"))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
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
