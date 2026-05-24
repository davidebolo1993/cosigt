rule impg_refine:
	'''
	https://github.com/pangenome/impg
	- Impg refine
	'''
	input:
		paf=get_merged_paf,
		bed=region_bed_path,
		flagger=rules.write_flagger_blacklist.output,
		index=rules.impg_index.output
	output:
		haplotypes_bed=outpath("refine/impg/{chr}/{region}/{region}.haplotypes.bed"),
		refined_bed=outpath("refine/impg/{chr}/{region}/{region}.refined.bed")
	threads:
		4
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['mid']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['mid']['runtime']
	container:
		'docker://davidebolo1993/impg:0.3.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.impg_refine.benchmark.txt'
	params:
		pansn=config['pansn_prefix']
	shell:
		'''
		impg \
			refine \
			-p {input.paf} \
			-b <(awk -v var={params.pansn} 'NF>=4{{print var$1,$2,$3,$4}} NF==3{{print var$1,$2,$3}}' OFS="\\t" {input.bed}) \
			-i {input.index} \
			-d 200000 \
			--span-bp 1000 \
			--pansn-mode haplotype \
			--extension-step 10000 \
			--support-output {output.haplotypes_bed} \
			--blacklist-bed {input.flagger} \
			-t 4 \
			> {output.refined_bed}
		'''

def get_all_optimal_beds(wildcards):
	'''
	https://github.com/davidebolo1993/cosigt
	- Generate all optimal bed files from impg_refine rule
	'''
	return [
		outpath("refine", "impg", REGION_ROWS[region]["chrom"], region, f"{region}.refined.bed")
		for region in REGION_ORDER
	]

rule make_bed_refined:
	'''
	https://github.com/davidebolo1993/cosigt
	- Output the final bed with all the regions refined
	'''
	input:
		get_all_optimal_beds
	output:
		outpath("refine/regions_refined.bed")
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default']['small']['mem_mb'],
		runtime=lambda wildcards, attempt: attempt * config['default']['small']['runtime']
	benchmark:
		'benchmarks/make_bed_refined.benchmark.txt'
	shell:
		'''
		for file in {input}; do
			awk 'NR>1 {{split($1, a, "#"); print a[3], $2, $3, $4}}' OFS="\\t" "$file"
		done | sort -k1,1 -k2,2n | uniq > {output}
		'''
