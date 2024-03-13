from glob import glob
import os

rule odgi_build:
	'''
	odgi build
	'''
	input:
		config['graph']
	output:
		config['output'] + '/odgi/build/' + os.path.basename(config['graph'].replace('.gfa', '.og'))
	threads:
		config['odgi']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['odgi']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['odgi']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	shell:
		'''
		odgi build \
		-g {input} \
		-o {output} \
		-t {threads}
		'''	

rule odgi_extract:
	'''
	odgi extract
	'''
	input:
		graph=rules.odgi_build.output,
		region=lambda wildcards: glob('resources/regions/{region}.bed'.format(region=wildcards.region))
	output:
		config['output'] + '/odgi/extract/{region}.og'
	threads:
		config['odgi']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['odgi']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['odgi']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	shell:
		'''
		odgi extract \
		-i {input.graph} \
		-b {input.region} \
		-t {threads} \
		-O \
		-o {output}
		'''


rule odgi_paths_fasta:
	'''
	odgi paths
	'''
	input:
		rules.odgi_extract.output
	output:
		config['output'] + '/odgi/paths/fasta/{region}.fa'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://pangenome/odgi:1707818641'
	shell:
		'''
		odgi paths \
		-i {input} \
		-f > {output}
		'''

rule faidx:
	'''
	samtools faidx
	'''
	input:
		rules.odgi_paths_fasta.output
	output:
		config['output'] + '/odgi/paths/fasta/{region}.fa.fai'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	shell:
		'''
		samtools faidx {input}
		'''


rule get_proper_region:
	input:
		rules.faidx.output
	output:
		config['output'] + '/samtools/bed/{region}.bed'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default']['time']
	params:
		path=config['path']
	shell:
		'''
		grep {params.path} {input} | cut -f 1 | cut -d "#" -f 2 | tr ':' '\t'| tr '-' '\t' > {output}
		'''
	
rule samtools_view:
	'''
	samtools view to extract the region
	'''
	input:
		sample=lambda wildcards: glob('resources/alignment/{sample}.*am'.format(sample=wildcards.sample)),
		bed=rules.get_proper_region.output
	output:
		config['output'] + '/samtools/view/{sample}/{region}.bam'
	threads:
		config['samtools']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	params:
		ref=config['reference'],
		region='{region}'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	shell:
		'''
		samtools view \
		-O bam \
		-o {output} \
		-T {params.ref} \
		-@ {threads} \
		-L {input.bed} \
		-M \
		{input.sample}
		'''

rule samtools_fasta:
	'''
	samtools fasta to get fasta files from .bam
	'''
	input:
		rules.samtools_view.output
	output:
		config['output'] + '/samtools/fasta/{sample}/{region}.fasta.gz'
	threads:
		config['samtools']['threads']
	container:
		'docker://davidebolo1993/graph_genotyper:latest'
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['samtools']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['samtools']['time']
	shell:
		'''
		samtools fasta \
		-@ {threads} \
		{input} | pigz > {output}
		'''