# Workflow

Workflow for genotyping from graph

```bash
singularity pull docker://davidebolo1993/graph_genotyper:latest
singularity run graph_genotyper_latest.sif genotype.sh --help
```

# evaluation

```bash
./cosigt_idx-bwa-mem.sh /lizardfs/erikg/amylase_diversity_project/pggb/amy.41/*.gfa amy.41 16 /scratch
cat samples_todo.txt | while read s; do sbatch -c 4 --wrap './cosigt-bwa-gfainject.sh amy.41 /lizardfs/erikg/amylase/1kgseqs/'$s'-AMY.fa.gz amy.41_1kg.bwa-gfainject.6a5a145/'$s' 4'; done >x.jobids
cat amy.29_1kg.bwa-gfainject.6a5a145/*/*best* | grep -v ^# | tr -d '"' | strswap -i /dev/stdin -s hap.swaps| sed s/.gaf// | strswap -i /dev/stdin -s sample_hap_gts.txt | grep '#' | cut -f 2- -d'/' | tr _ '\t'  | tr '#' '\t' | tr '/' '\t' | awk '{ print $1; print $2; print $3; print $4; print $5; print " ";  sum+=(($2 == $4 || $2 == $5) && ($3 == $4 || $3 == $5)) } END { print sum, NR, sum/NR; }'
```