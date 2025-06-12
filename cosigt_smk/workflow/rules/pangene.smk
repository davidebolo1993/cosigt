rule pangene_getaa:
	'''
	https://github.com/lh3/pangene
	'''
	input:
		gtf=rules.subset_gtf.output,
		proteins=config['proteins']
	output:
		config['output'] + '/pangene/proteins/{chr}/{region}/proteins.faa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/pangene:1.1'
	conda:
		'../envs/pangene.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_getaa.benchmark.txt'
	shell:
		'''
		pangene.js \
			getaa \
			-c \
			{input.gtf} \
			{input.proteins} > {output}
		'''

checkpoint pangene_prepare:
	'''
	https://github.com/lh3/miniprot
	'''
	input:
		asm=rules.bedtools_getfasta.output,
		proteins=rules.pangene_getaa.output
	output:
		directory(config['output'] + '/pangene/assemblies/{chr}/{region}_single_assemblies')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/pangene:1.1'
	conda:
		'../envs/pangene.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_prepare.benchmark.txt'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		bash workflow/scripts/pangene_prepare.sh {input.asm} {input.proteins} {params.pansn} {output}
		'''

def get_subpafs(wildcards):
	'''
	extract sub-pafs
	'''
	chromosome=wildcards.chr
	region=wildcards.region
	checkpoint_output = checkpoints.pangene_prepare.get(chr=chromosome, region=region).output[0]
	subpafs_files = glob(checkpoint_output + '/*_oriented.paf')
	return [os.path.basename(f).split('_oriented')[0] for f in subpafs_files]

rule pangene_graph:
	'''
	https://github.com/lh3/pangene
	'''
	input:
		lambda wildcards: expand(config['output'] + '/pangene/assemblies/{chr}/{region}_single_assemblies/{asm}_oriented.paf', chr=wildcards.chr, region=wildcards.region, asm=get_subpafs(wildcards))
	output:
		config['output'] + '/pangene/assemblies/{chr}/{region}/pangene_graph.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/pangene:1.1'
	conda:
		'../envs/pangene.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_graph.benchmark.txt'
	shell:
		'''
		pangene {input} --bed > {output}
		'''

rule pangene_graph_to_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		rules.pangene_graph.output
	output:
		config['output'] + '/pangene/viz/{chr}/{region}/plot.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.pangene_graph_to_bed.benchmark.txt'
	shell:
		'''
		sh workflow/scripts/convert_bed.sh {input} > {output}
		'''

rule pangene_viz:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		bed=rules.pangene_graph_to_bed.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/pangene/viz/{chr}/{region}/genes.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_viz.benchmark.txt'
	shell:
		'''
		Rscript workflow/scripts/plotgggenes.r {input.bed} {input.json} {output}
		'''