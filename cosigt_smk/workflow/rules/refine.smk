rule impg_refine:
	'''
	https://github.com/pangenome/impg
	- Impg refine
	'''
	input:
		paf=get_merged_paf,
		bed=lambda wildcards: glob('resources/regions/{chr}/{region}.bed'.format(chr=wildcards.chr, region=wildcards.region))
	output:
		haplotypes_bed=config['output'] + '/refine/impg/{chr}/{region}/{region}.haplotypes.bed',
		refined_bed=config['output'] + '/refine/impg/{chr}/{region}/{region}.refined.bed'
	threads:
		4
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_mid']['time']
	container:
		'docker://davidebolo1993/impg:0.3.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.impg_refine.benchmark.txt'
	params:
		flagger_blacklist=config['flagger_blacklist'],
		pansn=config['pansn_prefix']
	shell:
		'''
		impg \
			refine \
			-p {input.paf} \
			-b <(awk -v var={params.pansn} 'NF>=4{{print var$1,$2,$3,$4}} NF==3{{print var$1,$2,$3}}' OFS="\\t" {input.bed}) \
			-d 200000 \
			--span-bp 1000 \
			--pansn-mode haplotype \
			--extension-step 10000 \
			--support-output {output.haplotypes_bed} \
			--blacklist-bed {params.flagger_blacklist} \
			-t 4 \
			> {output.refined_bed}
		'''

def get_all_optimal_beds(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- Generate all optimal bed files from impg_refine rule
	'''
	all_files = []
	with open(config['all_regions']) as f:
		for line in f:
			fields = line.rstrip().split('\t')
			chr_name = fields[0]
			start = fields[1]
			end = fields[2]
			region = '_'.join([chr_name, start, end])
			all_files.append(f"{config['output']}/refine/impg/{chr_name}/{region}/{region}.refined.bed")
	return all_files

rule make_bed_refined:
	'''
	https://github.com/davidebolo1993/cosigt
	- Output the final bed with all the regions refined
	'''
	input:
		get_all_optimal_beds
	output:
		config['output'] + '/refine/regions_refined.bed'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/make_bed_refined.benchmark.txt'
	shell:
		'''
		for file in {input}; do
			awk 'NR>1 {{split($1, a, "#"); print a[3], $2, $3, $4}}' OFS="\\t" "$file"
		done | sort -k1,1 -k2,2n | uniq > {output}
		'''
	