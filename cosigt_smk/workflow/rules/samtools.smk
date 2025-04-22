rule samtools_fasta_mapped:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		sample=lambda wildcards: glob('resources/alignments/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.make_reference_bed.output,
		fasta=config['reference']
	output:
		config['output'] + '/samtools/fasta/{sample}/{chr}/{region}.mapped.fasta.gz'
	threads:
		config['samtools']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'	
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.samtools_fasta_mapped.benchmark.txt'
	shell:
		'''
		samtools view \
		-T {input.fasta} \
		-@ {threads} \
		-L {input.bed} \
		-M \
		-b \
		{input.sample} | \
		samtools sort -n | \
		samtools fasta \
		-@ {threads} \
		- | gzip > {output}
		'''

rule samtools_faidx_extract:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		asm_fai=lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
		ref_fasta=rules.pansnspec_target.output,
		ref_fai=rules.samtools_faidx_target.output,
		asm_bed=rules.filter_outliers.output,
		ref_bed=rules.make_reference_bed.output
	output:
		config['output'] + '/samtools/faidx/{chr}/{region}.fasta'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	benchmark:
		'benchmarks/{chr}.{region}.samtools_faidx_extract.benchmark.txt'
	conda:
		'../envs/samtools.yaml'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		asm_fai=$(echo {input.asm_fai})
		asm_fasta=$(echo "${{asm_fai%.*}}")
		samtools faidx \
		-r <(awk '{{print $1":"$2"-"$3}}' {input.asm_bed}) \
		$asm_fasta > {output}
		samtools faidx \
		-r <(awk -v var={params.pansn} '{{print var$1":"$2"-"$3}}' {input.ref_bed}) \
		{input.ref_fasta} >> {output}
		'''

rule samtools_faidx_index:
	'''
	https://github.com/samtools/samtools
	'''
	input:
		rules.samtools_faidx_extract.output
	output:
		config['output'] + '/samtools/faidx/{chr}/{region}.fasta.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/samtools:1.21'
	conda:
		'../envs/samtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.samtools_faidx_index.benchmark.txt'
	shell:
		'''
		samtools faidx {input}
		'''