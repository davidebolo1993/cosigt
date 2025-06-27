#SHOULD BE ONLY USED FOR REGION REFINEMENT
rule generate_expanded_regions:
	'''
	https://github.com/davidebolo1993/cosigt
	- Expand given region on the right and on the left flank
	'''
	output:
		temp(expand(config['output'] + '/refine/regions_expanded/{{chr}}/{{region}}/flank{flank}.bed', flank=range(0, 60000, 10000)))
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

rule impg_project_batches_expanded:
	'''
	https://github.com/pangenome/impg
	- Same stuff as impg.smk but expanding for flanks
	'''
	input:
		paf=get_merged_paf,
		bed=lambda wildcards: config['output'] + f'/refine/regions_expanded/{wildcards.chr}/{wildcards.region}/flank{wildcards.flank}.bed'
	output:
		unfiltered=temp(config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/unfiltered.bedpe.gz'),
		noblck=temp(config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.bedpe.gz'),
		merged=temp(config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.merged.bedpe.gz'),
		filtered=temp(config['output'] + '/refine/impg/{chr}/{region}/flank{flank}/noblck.merged.filtered.bedpe.gz')
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_mid']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	container:
		'docker://davidebolo1993/impg:0.2.4'
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
			-b <(awk -v var={params.pansn} '{{print var$1,$2,$3}}' OFS="\\t" {input.bed}) | gzip > {output.unfiltered}
		(zgrep -v \
			-E \
			-f {params.blacklist} {output.unfiltered} | bedtools \
			intersect \
			-a - \
			-b {params.flagger_blacklist} \
			-F 0.1 \
			-v \
			-wa | bedtools sort -i - || true) | gzip > {output.noblck}
		(zcat {output.noblck} | sh workflow/scripts/bedpe_merge.sh - 200000 || true) | gzip > {output.merged}
		(zcat {output.merged} | sh workflow/scripts/bedpe_filter.sh - {params.region} 1000 || true) | gzip > {output.filtered}
		'''

rule count_haplotypes_per_flank:
	'''
	https://github.com/davidebolo1993/cosigt
	- Number of contigs survived
	'''
	input:
		rules.impg_project_batches_expanded.output.filtered
	output:
		temp(config['output'] + '/refine/haplotype_counts/{chr}/{region}/flank{flank}.count')
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt * config['default_small']['mem_mb'],
		time=lambda wildcards, attempt: attempt * config['default_small']['time']
	benchmark:
		'benchmarks/{chr}.{region}.{flank}.count_haplotypes_per_flank.txt'
	shell:
		'''
		if [ -s {input} ]; then
			zcat {input} | cut -f1 | sort -u | wc -l > {output}
		else
			echo "0" > {output}
		fi
		'''

rule find_optimal_flank:
	'''
	https://github.com/davidebolo1993/cosigt
	- Just output the optimal bed file for that region
	'''
	input:
		counts=expand(config['output'] + '/refine/haplotype_counts/{{chr}}/{{region}}/flank{flank}.count', flank=range(0, 60000, 10000)),
		beds=expand(config['output'] + '/refine/regions_expanded/{{chr}}/{{region}}/flank{flank}.bed', flank=range(0, 60000, 10000))
	output:
		optimal_bed=temp(config['output'] + '/refine/{chr}/{region}/{region}.refined.bed'),
		summary=config['output'] + '/refine/{chr}/{region}/{region}.summary.txt'
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

def get_all_optimal_beds(wildcards):
    '''
	https://github.com/davidebolo1993/cosigt
    - Generate all optimal bed files from find_optimal_flank rule
    '''
    all_files = []
    with open(config['all_regions']) as f:
        for line in f:
            fields = line.rstrip().split('\t')
            chr_name = fields[0]
            start = fields[1]
            end = fields[2]
            region = '_'.join([chr_name, start, end])
            all_files.append(f"{config['output']}/refine/{chr_name}/{region}/{region}.refined.bed")
    return all_files


rule make_bed_refined:
	'''
	https://github.com/davidebolo1993/cosigt
	- Intersect the chosen regions with the original regions
	- This is mainly used to get back the annotations if they were provided initially
	- Ouput the final bed
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
	container:
		'docker://davidebolo1993/bedtools:2.31.0'
	benchmark:
		'benchmarks/make_bed_refined.benchmark.txt'
	conda:
		'../envs/bedtools.yaml'
	params:
		regions=config['all_regions']
	shell:
		'''
		bedtools intersect \
			-a <(cat {input} | bedtools sort -i -)  \
			-b <(bedtools sort -i {params.regions}) \
			-wao | \
			cut -f 1,2,3,7 > {output}
		'''
	