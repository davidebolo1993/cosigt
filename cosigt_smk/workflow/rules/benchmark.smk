rule make_tpr_table:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		samples=lambda wildcards: expand(config['output'] + '/cosigt/{sample}/{chr}/{region}/sorted_combos.tsv', sample=config['samples'], chr='{chr}', region='{region}'),
		tsv=config['output'] + '/odgi/dissimilarity/{chr}/{region}.tsv',
		json=config['output'] + '/cluster/{chr}/{region}.clusters.json'
	output:
		config['output'] + '/benchmark/{chr}/{region}/tpr.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	params:
		tsv=config['output'] + '/cluster/{chr}/{region}.clusters.hapdist.tsv'
	shell:
		'''
		Rscript \
		workflow/scripts/calc_tpr.r \
		{input.tsv} \
		{input.json} \
		{params.tsv} \
		{output} \
		{input.samples}
		'''

rule odgi_flip_pggb_graph:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		og=rules.pggb_construct.output,
		bed=rules.make_reference_bed.output
		#bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		config['output'] + '/benchmark/{chr}/{region}/flip.og'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		odgi flip \
			-i {input.og} \
			-o {output} \
			--ref-flips <(awk -v var={params.pansn} '{{print var$1":"$2"-"$3}}' {input.bed})
		'''	

rule odgi_og_to_fasta:
	'''
	https://github.com/pangenome/odgi
	'''
	input:
		rules.odgi_flip_pggb_graph.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/flip.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://pangenome/odgi:1745375412'
	conda:
		'../envs/odgi.yaml'
	shell:
		'''
		odgi paths \
			-i {input} \
			-f | sed 's/_inv$//g'> {output}
		'''

checkpoint prepare_combinations_for_qv:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		tsv=rules.make_tpr_table.output,
		fasta=rules.odgi_og_to_fasta.output
	output:
		directory(config['output'] + '/benchmark/{chr}/{region}/qv_prep')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	shell:
		'''
		bash workflow/scripts/prepare_qv.sh {input.tsv} {input.fasta} {output}
		'''

def get_samples(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	import os
	chr=wildcards.chr
	region=wildcards.region
	checkpoint_output = checkpoints.prepare_combinations_for_qv.get(chr=chr, region=region).output[0]
	samples = [x[0] for x in os.walk(checkpoint_output)][1:]
	return [os.path.basename(f)for f in samples]

rule calculate_qv:
	'''
	https://github.com/Martinsos/edlib
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		config['output'] + '/benchmark/{chr}/{region}/qv_prep/{sample}/ids.tsv'
	output:
		config['output'] + '/benchmark/{chr}/{region}/qv_prep/{sample}/qv.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/edlib:1.2.7'
	params:
		outdir=config['output'] + '/benchmark/{chr}/{region}/qv_prep/{sample}'
	shell:
		'''
		bash workflow/scripts/calculate_qv.sh {params.outdir} {output}
		'''

rule combine_qv:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/benchmark/{chr}/{region}/qv_prep/{sample}/qv.tsv', chr='{chr}', region='{region}', sample=get_samples(wildcards))
	output:
		config['output'] + '/benchmark/{chr}/{region}/bestqv.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	shell:
		'''
		cat {input} > {output}
		'''

rule combine_tpr_qv:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		tpr=rules.make_tpr_table.output,
		qv=rules.combine_qv.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/tpr_qv.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	shell:
		'''
		Rscript \
		workflow/scripts/combine_tpr_qv.r \
		{input.tpr} \
		{input.qv} \
		{wildcards.region} \
		{output}
		'''

def get_all_tpr_qv_files(wildcards):
	'''
	all the files from the prev analysis
	'''
	all_files = []
	with open(config['all_regions']) as f:
		for line in f:
			fields = line.rstrip().split('\t')
			chr = fields[0]
			start= fields[1]
			end= fields[2]
			region='_'.join([chr, start, end])
			all_files.append(f"{config['output']}/benchmark/{chr}/{region}/tpr_qv.tsv")
	return all_files

rule plot_tpr:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		get_all_tpr_qv_files
	output:
		config['output'] + '/benchmark/tpr.edr.png',
		config['output'] + '/benchmark/tpr.qv.png',
		config['output'] + '/benchmark/tpr.tpr_bar.png',
		config['output'] + '/benchmark/tpr.qv_bar.png'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	params:
		annot_bed=config['all_regions'],
		output_prefix=config['output'] + '/benchmark/tpr'
	shell:
		'''
		Rscript \
		workflow/scripts/plot_tpr.r \
		{params.output_prefix} \
		{params.annot_bed} \
		{input}
		'''
