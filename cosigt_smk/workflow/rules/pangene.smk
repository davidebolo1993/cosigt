
rule subset_gtf:
	'''
	https://github.com/davidebolo1993/cosigt
	- Subset the gtf file to the region of interest
	'''
	input:
		gtf=config['gtf'],
		bed=rules.make_reference_bed.output
	output:
		config['output'] + '/annotations/{chr}/{region}/{region}.gtf.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['small']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['small']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.subset_gtf.benchmark.txt'
	shell:
		'''
		bedtools intersect \
		-wa \
		-a {input.gtf} \
		-b {input.bed} | gzip > {output}
		'''

rule pangene_getaa:
	'''
	https://github.com/lh3/pangene
	- Extract protein sequences from annotation
	'''
	input:
		gtf=rules.subset_gtf.output,
		proteins=config['proteins']
	output:
		config['output'] + '/pangene/proteins/{chr}/{region}/{region}.faa.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['mid']['time']
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
			{input.proteins} | bgzip > {output}
		'''

checkpoint pangene_prepare:
	'''
	https://github.com/lh3/miniprot
	- Protein-to-assemblies alignment using miniprot
	- Flip paf files if they do not match the reference cds orientation
	'''
	input:
		asm=rules.bedtools_getfasta.output.fasta,
		proteins=rules.pangene_getaa.output
	output:
		temp(directory(config['output'] + '/pangene/assemblies/{chr}/{region}/single_assemblies'))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['mid']['time']
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
	https://github.com/davidebolo1993/cosigt
	- Needed to get the wildcard name
	'''
	chromosome=wildcards.chr
	region=wildcards.region
	checkpoint_output = checkpoints.pangene_prepare.get(chr=chromosome, region=region).output[0]
	subpafs_files = glob(checkpoint_output + '/*_oriented.paf.gz')
	return [os.path.basename(f).split('_oriented')[0] for f in subpafs_files]

rule pangene_graph:
	'''
	https://github.com/lh3/pangene
	- Construct the pangene graph and convert to a suitable format
	'''
	input:
		lambda wildcards: expand(config['output'] + '/pangene/assemblies/{chr}/{region}/single_assemblies/{asm}_oriented.paf.gz', chr=wildcards.chr, region=wildcards.region, asm=get_subpafs(wildcards))
	output:
		config['output'] + '/pangene/assemblies/{chr}/{region}/{region}.plot.bed.gz'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt *  config['default']['high']['time']
	container:
		'docker://davidebolo1993/pangene:1.1'
	conda:
		'../envs/pangene.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_graph.benchmark.txt'
	shell:
		'''
		pangene {input} --bed | sh workflow/scripts/convert_bed.sh - | gzip > {output}
		'''

rule pangene_viz:
	'''
	https://github.com/davidebolo1993/cosigt
	- Plot gggenes-like visualisation of pangene results
	'''
	input:
		bed=rules.pangene_graph.output,
		fai=rules.bedtools_getfasta.output.fai,
		json=rules.make_clusters.output
	output:
		config['output'] + '/pangene/viz/{chr}/{region}/{region}.genes.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *  config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['high']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.pangene_viz.benchmark.txt'
	params:
		paf_folder=config['output'] + '/pangene/assemblies/{chr}/{region}/single_assemblies',
		tsv=config['output'] + '/cluster/{chr}/{region}/{region}.clusters.medoids.tsv'
	shell:
		'''
		Rscript workflow/scripts/plotgggenes.r {input.bed} {input.json} {input.fai} {output} {params.tsv}
		rm -rf {params.paf_folder}
		'''