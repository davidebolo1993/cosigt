rule gfa_inject:
	'''
	gfainject
	'''
	input:
		gfa=rules.odgi_view.output,
		bam=rules.bwa_mem2_samtools_sort.output
	output:
		config['output'] + '/gfainject/{sample}/{region}.gaf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/cosigt_workflow:latest'
	params:

	shell:
		'''
		gfainject --gfa {input.gfa} --bam {input.bam} > {output}
		'''