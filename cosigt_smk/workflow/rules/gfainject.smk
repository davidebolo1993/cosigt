rule gfa_inject:
	'''
	gfa_inject
	'''
	input:
		gfa=rules.odgi_chop.output,
		bam=rules.bwa_mem_samtools_sort.output
	output:
		"results/cosigt_results/{sample}/{sample}.x.gaf"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:

	shell:
		'''
		gfainject --gfa {input.gfa} --bam {input.bam} > {output}
		'''
