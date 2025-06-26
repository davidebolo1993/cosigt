rule make_reference_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	- Extract reference region based on bedpe output
	'''
	input:
		rules.impg_project_batches.output
	output:
		temp(config['output'] + '/samtools/bed/{chr}/{region}/{region}.bed.gz')
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_reference_bed.benchmark.txt'
	shell:
		'''
		zcat {input} | cut -f 4-6 | \
		bedtools sort -i - | \
		bedtools merge -i - | \
		sed 's/#/\t/g' | \
		rev | \
		cut -f 1-3 | \
		rev | gzip > {output}
		'''

rule bedtools_getfasta:
	'''
	https://github.com/arq5x/bedtools2
	- Extract the region of interest from the contigs
	- Do the same for the reference
	- Build index
	'''
	input:
		asm_fai=lambda wildcards: glob('resources/assemblies/{chr}/*fai'.format(chr=wildcards.chr)),
		ref_fasta=rules.pansnspec_target.output.fasta,
		ref_fai=rules.pansnspec_target.output.fai,
		asm_bed=rules.impg_project_batches.output.filtered,
		ref_bed=rules.make_reference_bed.output
	output:
		fasta=config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz',
		fai=config['output'] + '/bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fai'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_high']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	benchmark:
		'benchmarks/{chr}.{region}.bedtools_getfasta.benchmark.txt'
	conda:
		'../envs/bedtools.yaml'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		asm_fai=$(echo {input.asm_fai})
		asm_fasta=$(echo "${{asm_fai%.*}}")
		bedtools getfasta \
			-fi $asm_fasta \
			-bed {input.asm_bed} | bgzip -c > {output.fasta}
		bedtools getfasta \
			-fi {input.ref_fasta} \
			-bed <(zcat {input.ref_bed} | awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t") | bgzip -c >> {output.fasta}
		samtools faidx {output.fasta}
		'''
