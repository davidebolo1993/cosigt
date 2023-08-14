from glob import glob

rule get_proper_region:
	input:
		rules.faidx.output
	output:
		'results/samtools/bed/{region}.bed'
	threads:
		1
	params:
		path=config['path']
	shell:
		'''
		grep {params.path} {input} | cut -f 1 | cut -d "_" -f 1,3-4 | tr '_' '\t' > {output}
		'''
	
rule samtools_view:
	'''
	samtools view to extract the region
	'''
	input:
		sample=lambda wildcards: glob('resources/alignment/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.get_proper_region.output
	output:
		'results/samtools/view/{sample}/{region}.bam'
	threads:
		config['samtools']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		ref=config['reference'],
		region='{region}'
	resources:
		mem_mb=config['samtools']['mem_mb'],
		time=config['samtools']['time']
	shell:
		'''
		samtools view \
		-O bam \
		-o {output} \
		-T {params.ref} \
		-@ {threads} \
		-L {input.bed} \
		-M \
		{input.sample} \
		'''

rule samtools_fasta:
	'''
	samtools fasta to get fasta files from .bam
	'''
	input:
		rules.samtools_view.output
	output:
		'results/samtools/fasta/{sample}/{region}.fasta.gz'
	threads:
		config['samtools']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=config['samtools']['mem_mb'],
		time=config['samtools']['time']
	shell:
		'''
		samtools fasta \
		-@ {threads} \
		{input} | pigz > {output}
		'''
