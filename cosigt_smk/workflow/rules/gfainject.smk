rule gfainject_inject:
	'''
	https://github.com/ekg/gfainject
	'''
	input:
		gfa=rules.odgi_view.output,
		bam=rules.bwamem2_mem_samtools_sort.output
	output:
		config['output'] + '/gfainject/{sample}/{chr}/{region}.gaf.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/gfainject:0.2.0'
	conda:
		'../envs/gfainject.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.gfainject_inject.benchmark.txt'
	shell:
		'''
		gfainject \
		--gfa {input.gfa} \
		--bam {input.bam} \
		--alt-hits 10000 | gzip > {output}
		'''
