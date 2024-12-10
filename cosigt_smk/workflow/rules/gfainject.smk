rule gfainject_inject:
	'''
	https://github.com/ekg/gfainject
	'''
	input:
		gfa=rules.odgi_view.output,
		bam=rules.bwamem2_mem_samtools_sort.output
	output:
		config['output'] + '/gfainject/{sample}/{region}.gaf.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	conda:
		'../envs/gfainject.yaml'
	benchmark:
		'benchmarks/{sample}.{region}.gfainject_inject.benchmark.txt'
	shell:
		'''
		gfainject \
		--gfa {input.gfa} \
		--bam {input.bam} | gzip > {output}
		'''