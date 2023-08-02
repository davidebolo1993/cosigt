import pandas as pd

df=(pd.read_table(config['samples'], dtype={'sample_id': str, 'cram':str})
	.set_index('sample_id', drop=False)
	.sort_index()
)

rule cosigt_genotype:
	'''
	cosigt genotype
	'''
	input:
		zpath=rules.odgi_paths_matrix.output,
		xpack=rules.gafpack_coverage.output
	output:
		geno='results/cosigt/{sample}/{region}/best_genotype.tsv',
		combo='results/cosigt/{sample}/{region}/combos.tsv'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix='results/cosigt/{sample}/{region}'
	shell:
		'''
		cosigt {input.zpath} {input.xpack} resources/extra/bad_samples.txt {params.prefix}
		'''

rule evaluate_cosigt:
	'''
	cosigt evaluation
	'''
	input:
		rules.cosigt_genotype.output.combo
	output:
		'results/cosigt/{sample}/{region}/evaluation.tsv'
	threads:
		1
	params:
		samplename='{sample}'
	shell:
		'''
		sort -k 2 -n -r {input} | awk -v var="{params.samplename}" -F '{params.samplename}' '{{print NR "\t" NF-1 "\t" var "\t" $0}}' | sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | sed 's/[0-9]-[0-9]/_/g' |tr '-' '\t' | sed 's/haplotype1_/haplotype1-/g' | sed 's/haplotype2_/haplotype2-/g'  > {output}
		'''

rule plot_evaluation:
	'''
	plot evaluation. Need to wait for all the evaluations table to be there
	'''
	input:
		expand('results/cosigt/{sample}/{region}/evaluation.tsv', sample=df['sample_id'].tolist(),region=config['region']),
		lambda wildcards: glob('results/cosigt/*/{region}/evaluation.tsv'.format(region=wildcards.region))
	output:
		'results/cosigt/evaluation/{region}/evaluation.pdf'
	threads:
		1
	conda:
		'../envs/python.yml'
	params:
		region='{region}'
	shell:
		'''
		python workflow/scripts/tpr.py results/cosigt '' {output} {params.region}
		'''

rule get_best_n_clusters:
	'''
	cluster haplotypes by similarity
	'''
	input:
		mtx=rules.odgi_similarity.output
	output:
		'results/cosigt/evaluation/{region}/dendrogram.jaccard.bestcut.json'
	threads:
		1
	conda:
		'../envs/python.yml'
	params:
		region='{region}'	
	shell:
		'''
		python workflow/scripts/cluster.py {input.mtx} results/cosigt/evaluation/{params.region}
		'''

rule plot_evaluation_dendro_jaccard:
	'''
	plot evaluation using helper dendrogram - based on jaccard distance
	'''
	input:
		expand('results/cosigt/{sample}/{region}/evaluation.tsv', sample=df['sample_id'].tolist(),region=config['region']),
		lambda wildcards: glob('results/cosigt/*/{region}/evaluation.tsv'.format(region=wildcards.region)),
		json=rules.get_best_n_clusters.output
	output:
		'results/cosigt/evaluation/{region}/evaluation.jaccard.pdf'
	threads:
		1
	conda:
		'../envs/python.yml'
	params:
		region='{region}'
	shell:
		'''
		python workflow/scripts/tpr.py results/cosigt {input.json} {output} {params.region}
		'''