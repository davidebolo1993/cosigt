#SHOULD BE ONLY USED FOR REGION REFINEMENT
rule generate_expanded_regions:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	output:
		bed=expand(config['output'] + '/refine/regions_expanded/{{chr}}/{{region}}/flank{flank}.bed', flank=range(0, 60000, 10000))
	run:
		import os
		outdir=os.path.dirname(output[0])
		os.makedirs(outdir, exist_ok=True)
		base_chr, start, end = wildcards.region.split('_')
		start = int(start)
		end = int(end)
		for flank in range(0, 60000, 10000):
			new_start = max(0, start - flank)
			new_end = end + flank
			with open(f'{outdir}/flank{flank}.bed', 'w') as f:
				f.write(f'{base_chr}\t{new_start}\t{new_end}\n')

rule concatenate_paf_batches_per_region:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		lambda wildcards: expand(config['output'] + '/wfmash/{chr}/batches/paf/{batch}.paf', chr='{chr}', batch=get_batches(wildcards))
	output:
		config['output'] + '/wfmash/{chr}/merged/all_batches.paf'
	threads:
		1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_high']['time']
	benchmark:
		'benchmarks/{chr}.concatenate_paf_batches_per_region.txt'
	shell:
		'''
		cat {input} > {output}
		'''

rule impg_project_batches_expanded:
	'''
	https://github.com/pangenome/impg
	'''
	input:
		paf=rules.concatenate_paf_batches_per_region.output,
		bed=lambda wildcards: config['output'] + f'/refine/regions_expanded/{wildcards.chr}/{wildcards.region}/flank{wildcards.flank}.bed'
	output:
		unfiltered=config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/unfiltered.bedpe',
		noblck=config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.bedpe',
		merged=config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.merged.bedpe',
		filtered=config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.merged.filtered.bedpe'
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/impg:0.2.3'
	conda:
		'../envs/impg.yaml'
	benchmark:
		'benchmarks/{chr}.{region}.{flank}.impg_project_batches_expanded.benchmark.txt'
	params:
		blacklist=config['blacklist'],
		flagger_blacklist=config['flagger_blacklist'],
		pansn=config['pansn_prefix'],
		region='{region}'
	shell:
		'''
		impg \
		query \
		-p {input.paf} \
		-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\t" {input.bed}) > {output.unfiltered}
		(grep -v -E -f {params.blacklist} {output.unfiltered} | \
		bedtools intersect -a - -b {params.flagger_blacklist} -F 0.20 -v -wa | \
		bedtools sort -i - || true) > {output.noblck}
		(sh workflow/scripts/bedpe_merge.sh {output.noblck} 200000 || true) > {output.merged}
		(sh workflow/scripts/bedpe_filter_refine.sh {output.merged} {input.bed} 1000 || true) > {output.filtered}
		'''

rule count_haplotypes_per_flank:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.merged.filtered.bedpe'
	output:
		config['output'] + '/refine/haplotype_counts/{chr}/{region}/flank{flank}.count'
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.{flank}.count_haplotypes_per_flank.txt'
	shell:
		'''
		if [ -s {input} ]; then
			cut -f1 {input} | sort -u | wc -l > {output}
		else
			echo "0" > {output}
		fi
		'''

rule find_optimal_flank:
	'''
	https://github.com/davidebolo1993/cosigt
	'''
	input:
		counts=expand(config['output'] + '/refine/haplotype_counts/{{chr}}/{{region}}/flank{flank}.count', flank=range(0, 60000, 10000)),
		beds=expand(config['output'] + '/refine/regions_expanded/{{chr}}/{{region}}/flank{flank}.bed', flank=range(0, 60000, 10000))
	output:
		optimal_bed=config['output'] + '/refine/{chr}/{region}/roi_refined.bed',
		summary=config['output'] + '/refine/{chr}/{region}/refinement_summary.txt'
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.find_optimal_flank.txt'
	run:
		import os
		
		flank_counts = {}
		flanks = list(range(0, 60000, 10000))
		
		for i, flank in enumerate(flanks):
			count_file = input.counts[i]
			with open(count_file, 'r') as f:
				count = int(f.read().strip())
			flank_counts[flank] = count
		
		max_count = max(flank_counts.values())
		optimal_flank = min([flank for flank, count in flank_counts.items() if count == max_count])		
		optimal_bed_source = f"{config['output']}/refine/regions_expanded/{wildcards.chr}/{wildcards.region}/flank{optimal_flank}.bed"
		os.makedirs(os.path.dirname(output.optimal_bed), exist_ok=True)
		with open(optimal_bed_source, 'r') as src, open(output.optimal_bed, 'w') as dst:
			dst.write(src.read())
		#write a small summary
		with open(output.summary, 'w') as f:
			f.write(f"Region: {wildcards.region}\n")
			f.write(f"Chromosome: {wildcards.chr}\n")
			f.write(f"Optimal flank size: {optimal_flank}\n")
			f.write(f"Maximum haplotype count: {max_count}\n")
			f.write(f"Flank size breakdown:\n")
			for flank in sorted(flank_counts.keys()):
				f.write(f"  Flank {flank}: {flank_counts[flank]} haplotypes\n")
