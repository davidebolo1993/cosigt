rule cosigt_genotype:
	'''
	https://github.com/davidebolo1993/cosigt
	- This is the actual genotyping step
	'''
	input:
		graph_cov_map=rules.odgi_utils.output.paths,
		sample_cov_map=rules.gafpack_coverage.output,
		json=rules.make_clusters.output,
		mask=rules.filter_nodes.output
	output:
		geno=outpath("cosigt/{sample}/{chr}/{region}/{region}.cosigt_genotype.tsv"),
		combos=outpath("cosigt/{sample}/{chr}/{region}/{region}.sorted_combos.tsv.gz")
	group:
		"genotype"
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://davidebolo1993/cosigt:0.2'
	conda:
		'../envs/cosigt.yaml'
	params:
		prefix=outpath("cosigt/{sample}/{chr}/{region}"),
		sample_id='{sample}',
		ploidy=lambda wildcards: region_ploidy(wildcards.region)
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
			--ploidy {params.ploidy} \
			-m {input.mask}
		'''
	
rule samtools_faidx_besthaps_fasta:
	'''
	https://github.com/samtools/samtools
	- Get the predicted contigs out
	'''
	input:
		geno=rules.cosigt_genotype.output.geno,
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai
	output:
		fasta=temp(outpath("cosigt/{sample}/{chr}/{region}/viz/{region}.haplotypes.fasta")),
		regions=temp(outpath("cosigt/{sample}/{chr}/{region}/viz/{region}.haplotypes.regions.txt"))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	container:
		'docker://davidebolo1993/samtools:1.23.1'
	conda:
		'../envs/samtools.yaml'
	params:
		pansn=config['pansn_prefix'] + '{chr}',
		ploidy=lambda wildcards: region_ploidy(wildcards.region)
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.samtools_faidx_besthaps_fasta.benchmark.txt'
	shell:
		'''
		if [ {params.ploidy} -gt 2 ]; then
			echo "SVbyEye haplotype plotting supports ploidy 1 or 2; region {wildcards.region} has ploidy {params.ploidy}." >&2
			exit 1
		fi
		bash workflow/scripts/select_haplotypes.sh \
			{input.fai} \
			{input.geno} \
			all \
			{params.pansn} > {output.regions}
		samtools faidx -r {output.regions} {input.fasta} > {output.fasta}
		'''

rule minimap2_ava:
	'''
	https://github.com/lh3/minimap2
	- Realign predicted haplotypes, all-vs-all alignment
	'''
	input:
		rules.samtools_faidx_besthaps_fasta.output.fasta
	output:
		temp(outpath("cosigt/{sample}/{chr}/{region}/viz/{region}.haplotypes.paf")),
	threads:
		config['minimap2']['ava']['threads']
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['ava']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['minimap2']['ava']['runtime']
	container:
		'docker://davidebolo1993/minimap2:2.31'
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
			-t {threads} \
			{input} \
			{input} > {output}
		'''

rule plot_ava:
	'''
	https://github.com/daewoooo/SVbyEye
	- Plot the relationship between alignments
	'''
	input:
		rules.minimap2_ava.output
	output:
		outpath("cosigt/{sample}/{chr}/{region}/viz/{region}.ava.png"),
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['high']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['high']['runtime']
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
		if [ -s {input} ]; then
			Rscript \
				workflow/scripts/plotava.r \
				{input} \
				{output} \
				{params.pansn}
		else
			touch {output}
		fi
		'''

#OPTIONALLY CALL SVs FROM COSIGT-PREDICTED HAPLOTYPES:

rule minimap2_align_sort_haps:
	'''
	https://github.com/lh3/minimap2
	https://github.com/samtools/samtools
	- Align cosigt-predicted haplotypes to the reference genome; sort and index BAMs for downstream svim-asm
	'''
	input:
		geno=rules.cosigt_genotype.output.geno,
		fasta=rules.bedtools_getfasta.output.fasta,
		fai=rules.bedtools_getfasta.output.fai,
		ref=config['reference']
	output:
		hap1_bam=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap1.sorted.bam")),
		hap1_csi=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap1.sorted.bam.csi")),
		hap1_regions=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap1.regions.txt")),
		hap2_bam=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap2.sorted.bam")),
		hap2_csi=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap2.sorted.bam.csi")),
		hap2_regions=temp(outpath("cosigt/{sample}/{chr}/{region}/svim_asm/{region}.hap2.regions.txt"))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['minimap2']['ava']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['minimap2']['ava']['runtime']
	container:
		'docker://davidebolo1993/minimap2:2.31'
	conda:
		'../envs/minimap2.yaml'
	params:
		ploidy=lambda wildcards: region_ploidy(wildcards.region)
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.minimap2_align_sort_haps.benchmark.txt'
	shell:
		'''
		if [ {params.ploidy} -gt 2 ]; then
			echo "svim-asm supports ploidy 1 or 2; region {wildcards.region} has ploidy {params.ploidy}." >&2
			exit 1
		fi
		bash workflow/scripts/select_haplotypes.sh \
			{input.fai} \
			{input.geno} \
			haplotype.1 > {output.hap1_regions}
		samtools faidx \
			-r {output.hap1_regions} \
			{input.fasta} | \
			minimap2 -a -x asm20 --cs -r2k -t {threads} {input.ref} - | \
			samtools sort -o {output.hap1_bam} --write-index

		if [ {params.ploidy} -eq 2 ]; then
			bash workflow/scripts/select_haplotypes.sh \
				{input.fai} \
				{input.geno} \
				haplotype.2 > {output.hap2_regions}
			samtools faidx \
				-r {output.hap2_regions} \
				{input.fasta} | \
				minimap2 -a -x asm20 --cs -r2k -t {threads} {input.ref} - | \
				samtools sort -o {output.hap2_bam} --write-index
		else
			: > {output.hap2_regions}
			touch {output.hap2_bam} {output.hap2_csi}
		fi
		'''

rule svim_asm:
	'''
	https://github.com/eldariont/svim-asm
	- Call structural variants from cosigt-predicted haplotypes aligned to reference
	'''
	input:
		hap1_bam=rules.minimap2_align_sort_haps.output.hap1_bam,
		hap1_csi=rules.minimap2_align_sort_haps.output.hap1_csi,
		hap2_bam=rules.minimap2_align_sort_haps.output.hap2_bam,
		hap2_csi=rules.minimap2_align_sort_haps.output.hap2_csi,
		ref=config['reference']
	output:
		vcf=outpath("cosigt/{sample}/{chr}/{region}/svim_asm/variants.vcf")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/svim-asm:1.0.3'
	conda:
		'../envs/svim-asm.yaml'
	params:
		workdir=outpath("cosigt/{sample}/{chr}/{region}/svim_asm"),
		ploidy=lambda wildcards: region_ploidy(wildcards.region)
	benchmark:
		'benchmarks/{sample}.{chr}.{region}.svim_asm.benchmark.txt'
	shell:
		'''
		if [ {params.ploidy} -eq 1 ]; then
			svim-asm haploid \
				{params.workdir} \
				{input.hap1_bam} \
				{input.ref}
		elif [ {params.ploidy} -eq 2 ]; then
			svim-asm diploid \
				{params.workdir} \
				{input.hap1_bam} \
				{input.hap2_bam} \
				{input.ref}
		else
			echo "svim-asm supports ploidy 1 or 2; region {wildcards.region} has ploidy {params.ploidy}." >&2
			exit 1
		fi
		'''

##OPTIONALLY, MAKE A VCF WITH COSIGT GENOTYPES FOR THIS REGION

rule make_region_vcf:
	'''
	https://github.com/davidebolo1993/cosigt
	- Collects cosigt genotype TSVs for all samples in a region and writes
	  a single per-region VCF. Allele 0 is the reference path; all other
	  haplotypes are numbered 1..N-1 in alphabetical order and listed in
	  INFO/ALLELES. Each sample gets a phased GT column with one allele
	  index per configured ploidy level.
	'''
	input:
		tsv=lambda wildcards: expand(
			outpath(
				"cosigt",
				"{sample}",
				"_".join(wildcards.region.split("_")[:-2]),
				"{region}",
				"{region}.cosigt_genotype.tsv",
			),
			sample=config['samples'],
			region=wildcards.region
		),
		fai=lambda wildcards: expand(
			rules.bedtools_getfasta.output.fai,
			chr='_'.join(wildcards.region.split('_')[:-2]),
			region=wildcards.region
		),
		ref_fai=config['reference'] + '.fai'
	output:
		vcf=temp(outpath("cosigt/vcf/{region}.vcf"))
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	conda:
		'../envs/python.yaml'
	container:
		'docker://davidebolo1993/pythonenv:3.13.3'
	params:
		pansn=lambda wildcards: config['pansn_prefix'] + '_'.join(wildcards.region.split('_')[:-2])
	benchmark:
		'benchmarks/{region}.make_region_vcf.benchmark.txt'
	shell:
		'''
		python workflow/scripts/make_region_vcf.py \
			--fai     {input.fai} \
			--tsv     {input.tsv} \
			--output  {output.vcf} \
			--pansn   {params.pansn} \
			--ref-fai {input.ref_fai}
		'''

rule merge_sort_vcf:
	'''
	https://github.com/samtools/bcftools
	- Concatenates all per-region VCFs (non-overlapping, same sample set),
	  sorts by coordinate, bgzips and tabix-indexes the result.
	'''
	input:
		expand(
			outpath("cosigt/vcf/{region}.vcf"),
			region=config['regions']
		)
	output:
		vcf=outpath("cosigt/vcf/cosigt.vcf.gz"),
		tbi=outpath("cosigt/vcf/cosigt.vcf.gz.tbi")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/bcftools:1.23.1'
	conda:
		'../envs/bcftools.yaml'
	benchmark:
		'benchmarks/merge_sort_vcf.benchmark.txt'
	shell:
		'''
		bcftools concat \
			{input} | \
		bcftools sort \
			-O z \
			-o {output.vcf} \
		&& tabix -p vcf {output.vcf}
		'''
