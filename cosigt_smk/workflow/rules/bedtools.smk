if ALLELE_SOURCE == 'assemblies':

	rule make_reference_bed:
		'''
		https://github.com/davidebolo1993/cosigt
		- Extract the reference region from impg .bedpe output
		'''
		input:
			rules.impg_project_batches.output
		output:
			outpath("bedtools/reference_bed/{chr}/{region}/{region}.bed.gz")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
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

else:

	rule make_reference_bed:
		'''
		https://github.com/davidebolo1993/cosigt
		- Get the reference BED directly from user-provided regions for custom alleles
		'''
		input:
			region=region_bed_path
		output:
			outpath("bedtools/reference_bed/{chr}/{region}/{region}.bed.gz")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
		benchmark:
			'benchmarks/{chr}.{region}.make_reference_bed.benchmark.txt'
		params:
			chrom='{chr}'
		shell:
			'''
			awk -v chrom={params.chrom:q} '$1 == chrom' {input.region} | gzip > {output}
			'''


rule make_alignment_bed:
	'''
	https://github.com/davidebolo1993/cosigt
	- Make alignment .bed for extracting reads from the original .cram/.bam files
	- This can be the reference .bed if no alternative intervals were provided
	'''
	input:
		ref_bed=rules.make_reference_bed.output,
		ori_bed=region_bed_path
	output:
		outpath("bedtools/alignment_bed/{chr}/{region}/{region}.bed.gz")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	conda:
		'../envs/bedtools.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.make_alignment_bed.benchmark.txt'
	shell:
		'''
		bedtools intersect \
			-a {input.ref_bed} \
			-b <(bedtools sort -i {input.ori_bed}) \
			-nonamecheck \
			-u | gzip > {output}
		bedtools intersect \
			-a  <(bedtools sort -i {input.ori_bed}) \
			-b {input.ref_bed} \
			-nonamecheck \
			-v | gzip >> {output}
		'''


if ALLELE_SOURCE == 'assemblies':

	rule bedtools_getfasta:
		'''
		https://github.com/arq5x/bedtools2
		- Extract allele sequences from impg .bedpe output and append the reference
		- Compress with bgzip
		- Build index
		'''
		input:
			asm_fasta=assembly_fasta_path,
			asm_fai=assembly_fai_path,
			ref_fasta=rules.pansnspec_target.output.fasta,
			ref_fai=rules.pansnspec_target.output.fai,
			asm_bed=rules.impg_project_batches.output.filtered,
			ref_bed=rules.make_reference_bed.output
		output:
			fasta=outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz"),
			fai=outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fai")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['default']['high']['runtime']
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
			bedtools getfasta \
				-fi {input.asm_fasta} \
				-bed {input.asm_bed} | bgzip -c > {output.fasta}
			bedtools getfasta \
				-fi {input.ref_fasta} \
				-bed <(zcat {input.ref_bed} | awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t") | bgzip -c >> {output.fasta}
			samtools faidx {output.fasta}
			'''

else:

	rule bedtools_getfasta:
		'''
		https://github.com/davidebolo1993/cosigt
		- Normalize a custom per-region allele FASTA to the standard workflow path
		- Compress with bgzip if needed
		- Build index
		'''
		input:
			custom=custom_allele_fasta_path
		output:
			fasta=outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz"),
			fai=outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.fai"),
			gzi=outpath("bedtools/getfasta/{chr}/{region}/{region}.fasta.gz.gzi")
		threads:
			1
		resources:
			mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
			runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
		container:
			'docker://davidebolo1993/samtools:1.22'
		benchmark:
			'benchmarks/{chr}.{region}.bedtools_getfasta.benchmark.txt'
		conda:
			'../envs/samtools.yaml'
		shell:
			'''
			if [[ {input.custom:q} == *.gz ]]; then
				cp {input.custom} {output.fasta}
			else
				bgzip -c {input.custom} > {output.fasta}
			fi
			samtools faidx {output.fasta}
			'''
