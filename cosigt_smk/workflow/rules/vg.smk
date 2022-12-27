rule vg_index:
	'''
	vg index
	'''
	input:
		rules.odgi_chop.output
	output:
		gbz="resources/vg/index.giraffe.gbz",
		min="resources/vg/index.min",
		dist="resources/vg/index.dist"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix="resources/vg/index",
		prefix_dir="resources/vg"
	shell:
		'''
		vg \
		autoindex \
		-w giraffe \
		-g {input} \
		-t {threads} \
		-p {params.prefix} \
		-T {params.prefix_dir}
		'''

rule vg_giraffe:
	'''
	vg giraffe
	'''
	input:
		gbz=rules.vg_index.output.gbz,
		min=rules.vg_index.output.min,
		dist=rules.vg_index.output.dist,
		fastq=rules.samtools_fastq.output
	output:
		"results/cosigt_results/{sample}/{sample}.x.gaf"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		vg \
		giraffe \
		-Z {input.gbz} \
		-m {input.min} \
		-d {input.dist} \
		-f {input.fastq} \
		-o gaf > {output}
		'''
