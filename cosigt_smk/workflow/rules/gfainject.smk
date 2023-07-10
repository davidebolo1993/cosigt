rule gfa_inject:
	'''
	gfainject
	'''
	input:
		gfa=rules.odgi_chop.output,
		bam=rules.bwa_mem2_samtools_sort.output
	output:
		'results/cosigt_results/{sample}/{sample}.x.gaf'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:

	shell:
		'''
		gfainject --gfa {input.gfa} --bam {input.bam} > {output}
		'''
