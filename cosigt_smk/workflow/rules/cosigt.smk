import pandas as pd

df=(pd.read_table(config['samples'], dtype={"sample_id": str, "cram":str})
	.set_index("sample_id", drop=False)
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
		geno="results/cosigt_results/{sample}/best_genotype.tsv",
		combo="results/cosigt_results/{sample}/combos.tsv"
	threads:
		1
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		prefix="results/cosigt_results/{sample}"
	shell:
		'''
		cosigt {input.zpath} {input.xpack} {params.prefix}
		'''

rule evaluate_cosigt:
	'''
	cosigt evaluation
	'''
	input:
		rules.cosigt_genotype.output.combo
	output:
		"results/cosigt_results/{sample}/evaluation.tsv"
	threads:
		1
	params:
		samplename="{sample}"
	shell:
		'''
		sort -k 2 -n -r {input} | awk -v var="{params.samplename}" -F '{params.samplename}' '{{print NR "\t" NF-1 "\t" var}}'  > {output}  		
		'''

rule plot_evaluation:
	'''
	plot evaluation
	'''
	input:
		expand("results/cosigt_results/{sample}/evaluation.tsv", sample=df['sample_id'].tolist())
	output:
		"results/cosigt_results/evaluation.pdf"
	threads:
		1
	conda:
		"../envs/r.yml"
	shell:
		'''
		Rscript workflow/scripts/plot_evaluation.r results/cosigt_results 10 {output}
		'''
