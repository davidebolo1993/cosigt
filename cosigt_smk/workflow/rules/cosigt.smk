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
		zpath=rules.odgi_build.output,
		xpack=rules.gafpack_coverage.output
	output:
		geno='results/cosigt_results/{sample}/best_genotype.tsv',
		combo='results/cosigt_results/{sample}/combos.tsv'
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix='results/cosigt_results/{sample}'
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
		'results/cosigt_results/{sample}/evaluation.tsv'
	threads:
		1
	params:
		samplename='{sample}'
	shell:
		'''
		sort -k 2 -n -r {input} | awk -v var="{params.samplename}" -F '{params.samplename}' '{{print NR "\t" NF-1 "\t" var "\t" $0}}' | sed 's/haplotype1-/haplotype1_/g' | sed 's/haplotype2-/haplotype2_/g' | tr '-' '\t' | sed 's/haplotype1_/haplotype1-/g' | sed 's/haplotype2_/haplotype2-/g'  > {output}
		'''

rule plot_evaluation:
	'''
	plot evaluation
	'''
	input:
		expand('results/cosigt_results/{sample}/evaluation.tsv', sample=df['sample_id'].tolist())
	output:
		'results/cosigt_results/evaluation/evaluation.pdf'
	threads:
		1
	conda:
		'../envs/python.yml'
	shell:
		'''
		python workflow/scripts/tpr.py results/cosigt_results '' {output}
		'''

rule get_best_n_clusters:
	'''
	cluster haplotypes by similarity
	'''
	input:
		files=expand('results/cosigt_results/{sample}/evaluation.tsv', sample=df['sample_id'].tolist()),
		mtx=rules.odgi_similarity.output
	output:
		'results/cosigt_results/evaluation/dendrogram.jaccard.bestcut.json'
	threads:
		1
	conda:
		'../envs/python.yml'
	shell:
		'''
		python workflow/scripts/cluster.py {input.mtx} results/cosigt_results/evaluation
		'''

rule plot_evaluation_dendro_jaccard:
	'''
	plot evaluation using helper dendrogram - based on jaccard distance
	'''
	input:
		rules.get_best_n_clusters.output
	output:
		'results/cosigt_results/evaluation/evaluation.jaccard.pdf'
	threads:
		1
	conda:
		'../envs/python.yml'
	shell:
		'''
		python workflow/scripts/tpr.py results/cosigt_results {input} {output}
		'''