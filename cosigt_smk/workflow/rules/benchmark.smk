rule make_tpr_table:
	'''
	https://github.com/davidebolo1993/cosigt
	- Compute the tpr
	- +1 TP if both the predicted haplotypes falls in the ground-truth haplotype clusters
	- +1 FN if at least one of the predicted haplotypes does not follow TP logic
	'''
	input:
		samples=lambda wildcards: expand(config['output'] + '/cosigt/{sample}/{chr}/{region}/{region}.sorted_combos.tsv.gz', sample=config['samples'], chr='{chr}', region='{region}'),
		tsv=rules.odgi_dissimilarity.output,
		json=rules.make_clusters.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/tpr.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['high']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	params:
		tsv=config['output'] + '/cluster/{chr}/{region}/{region}.clusters.hapdist.tsv'
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

rule odgi_flip_pggb_graph_to_fasta:
	'''
	https://github.com/pangenome/odgi
	- Orient the haplotypes with respect to the target
	- Output fasta file
	'''
	input:
		rules.pggb_construct.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/{region}.flipped.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['small']['time']
	container:
		'docker://pangenome/odgi:1753347183'
	conda:
		'../envs/odgi.yaml'
	params:
		pansn=config['pansn_prefix'],
		prefix=config['output'] + '/benchmark/{chr}/{region}'
	shell:
		'''
		odgi paths \
			-i {input} \
			-L | grep {params.pansn} > {params.prefix}/ref_path.txt
		odgi flip \
			-i {input} \
			-o - \
			--ref-flips <(cat {params.prefix}/ref_path.txt) | \
		odgi paths \
			-i - \
			-f | sed 's/_inv$//g' > {output}
		rm {params.prefix}/ref_path.txt
		'''	

rule prepare_combinations_for_qv:
	'''
	https://github.com/samtools/samtools
	- Prepare all possible combinations for QV calculation
	'''
	input:
		tsv=rules.make_tpr_table.output,
		fasta=rules.odgi_flip_pggb_graph_to_fasta.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/qv_prep.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['small']['time']
	container:
		'docker://davidebolo1993/samtools:1.22'
	conda:
		'../envs/samtools.yaml'
	params:
		outdir=config['output'] + '/benchmark/{chr}/{region}/qv'
	shell:
		'''
		if [ -f {output} ]; then
			rm {output}
		fi
		if [ -d {params.outdir} ]; then
			rm -rf {params.outdir}
		fi
		bash workflow/scripts/prepare_qv.sh {input.tsv} {input.fasta} {params.outdir} \
		&& touch {output}
		'''

rule calculate_qv:
	'''
	https://github.com/Martinsos/edlib
	https://github.com/davidebolo1993/cosigt
	- Actually, calculate the qv
	'''
	input:
		rules.prepare_combinations_for_qv.output
	output:
		config['output'] + '/benchmark/{chr}/{region}/qv_calc.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['high']['time']
	container:
		'docker://davidebolo1993/edlib:1.2.7'
	params:
		indir=config['output'] + '/benchmark/{chr}/{region}/qv'
	shell:
		'''
		if [ -f {output} ]; then
			rm {output}
		fi
		bash workflow/scripts/calculate_qv.sh {params.indir} \
		&& touch {output} \
		&& rm {params.indir}/*/qv.tmp.tsv
		'''

rule combine_tpr_qv:
	'''
	https://github.com/davidebolo1993/cosigt
	- Combine tpr and all qvs for futher plotting
	'''
	input:
		tpr=rules.make_tpr_table.output,
		qv=lambda wildcards: expand(config['output'] + '/benchmark/{chr}/{region}/qv_calc.done', chr='{chr}', region='{region}')
	output:
		config['output'] + '/benchmark/{chr}/{region}/{region}.tpr_qv.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['small']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	params:
		indir=config['output'] + '/benchmark/{chr}/{region}/qv',
		outqv=config['output'] + '/benchmark/{chr}/{region}/bestqv.tsv'
	shell:
		'''
		cat {params.indir}/*/qv.tsv > {params.outqv}
		Rscript \
			workflow/scripts/combine_tpr_qv.r \
			{input.tpr} \
			{params.outqv} \
			{wildcards.region} \
			{output}
		'''

def get_all_tpr_qv_files(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	Get all the files from the previuos analysis
	'''
	all_files = []
	with open(config['all_regions']) as f:
		for line in f:
			fields = line.rstrip().split('\t')
			chr = fields[0]
			start= fields[1]
			end= fields[2]
			region='_'.join([chr, start, end])
			all_files.append(f"{config['output']}/benchmark/{chr}/{region}/{region}.tpr_qv.tsv")
	return all_files

rule plot_tpr:
	'''
	https://github.com/davidebolo1993/cosigt
	- Make final plot
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
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb']
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
	