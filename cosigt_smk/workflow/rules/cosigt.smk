rule cosigt_genotype:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		graph_cov_map=rules.odgi_paths_matrix.output,
		sample_cov_map=rules.gafpack_coverage.output,
		json=rules.make_clusters.output,
		mask=rules.filter_nodes.output
	output:
		geno=config['output'] + '/cosigt/{sample}/{chr}/{region}/cosigt_genotype.tsv',
		combos=config['output'] + '/cosigt/{sample}/{chr}/{region}/sorted_combos.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/cosigt:0.1.3'
	conda:
		'../envs/cosigt.yaml'
	params:
		prefix=config['output'] + '/cosigt/{sample}/{chr}/{region}',
		sample_id='{sample}'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.cosigt_genotype.benchmark.txt'
	shell:
		'''
		cosigt \
		-p {input.graph_cov_map} \
		-g {input.sample_cov_map} \
		-c {input.json} \
		-o {params.prefix} \
		-i {params.sample_id} \
		-m {input.mask}
		'''

rule cosigt_genotype_submasks:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		graph_cov_map=rules.odgi_paths_matrix.output,
		sample_cov_map=rules.gafpack_coverage.output,
		json=rules.make_clusters.output,
		mask=config['output'] + '/odgi/paths/matrix/{chr}/{region}_submasks/{num}.tsv'
	output:
		geno=config['output'] + '/cosigt/{sample}/{chr}/{region}/masks/{num}/cosigt_genotype.tsv',
		combos=config['output'] + '/cosigt/{sample}/{chr}/{region}/masks/{num}/sorted_combos.tsv'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/cosigt:0.1.3'
	conda:
		'../envs/cosigt.yaml'
	params:
		prefix=config['output'] + '/cosigt/{sample}/{chr}/{region}/masks/{num}',
		sample_id='{sample}'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.{num}.cosigt_genotype_subregions.benchmark.txt'
	shell:
		'''
		cosigt \
		-p {input.graph_cov_map} \
		-g {input.sample_cov_map} \
		-c {input.json} \
		-o {params.prefix} \
		-i {params.sample_id} \
		-m {input.mask}
		'''

rule collect_submasks_genotypes:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/cosigt/{sample}/{chr}/{region}/masks/{num}/cosigt_genotype.tsv', chr=wildcards.chr, region=wildcards.region, sample=wildcards.sample, num=get_submasks(wildcards))
	output:
		config['output'] + '/cosigt/{sample}/{chr}/{region}/submasks.done'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.collect_submasks_genotypes.benchmark.txt'
	shell:
		'''
		touch {output}
		'''
	
rule samtools_faidx_besthaps_fasta:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		geno=rules.cosigt_genotype.output.geno,
		fasta=rules.samtools_faidx_extract.output,
		fai=rules.samtools_faidx_index.output
	output:
		config['output'] + '/cosigt/{sample}/{chr}/{region}/viz_files/haplotypes.fasta',
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	params:
		pansn=config['pansn_prefix'] + '{chr}'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.samtools_faidx_besthaps_fasta.benchmark.txt'
	shell:
		'''
		samtools faidx -r <(grep {params.path} {input.fai} | awk '{{print $1":1-"$2}}') {input.fasta} > {output} \
		&& samtools faidx -r <(grep $(cut -f 2 {input.geno} | tail -1) {input.fai} | awk '{{print $1":1-"$2}}') {input.fasta} >> {output} \
		&& samtools faidx -r <(grep $(cut -f 3 {input.geno} | tail -1) {input.fai} | awk '{{print $1":1-"$2}}') {input.fasta} >> {output}
		'''

rule minimap2_ava:
	'''
	https://github.com/lh3/minimap2
	'''
	input:
		rules.samtools_faidx_besthaps_fasta.output
	output:
		config['output'] + '/cosigt/{sample}/{chr}/{region}/viz_files/haplotypes.paf',
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['minimap2']['time']
	container:
		'docker://davidebolo1993/minimap2:2.28'
	conda:
		'../envs/minimap2.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.minimap2_allvall.benchmark.txt'
	shell:
		'''
		minimap2 \
			-x asm20 \
			-c \
			--eqx \
			-D \
			-P \
			--dual=no \
			{input} \
			{input} > {output}
		'''

rule plot_ava:
	'''
	https://github.com/daewoooo/SVbyEye
	'''
	input:
		rules.minimap2_ava.output
	output:
		config['output'] + '/cosigt/{sample}/{chr}/{region}/ava.pdf',
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/renv:4.3.3'
	conda:
		'../envs/r.yaml'
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.plot_ava.benchmark.txt'
	params:
		pansn=config['pansn_prefix'] + '{chr}'
	shell:
		'''
		Rscript \
			workflow/scripts/plotava.r \
			{input} \
			{output} \
			{params.path}
		'''