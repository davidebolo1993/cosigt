# cosigt

cosigt is a tool for genotyping complex loci in pangenome graphs. It takes as input a path coverage matrix for a locus from a pangenome graph (generated with odgi paths) and read alignments to the graph locus for a sample (converted to graph alignments with gfainject). It then compares all combinations of two haplotypes from the graph to the sample alignments to find the best matching pair of haplotypes.

## Usage

```
cosigt path_matrix.tsv sample.gafpack.gz bad_samples.txt prefix
```

- path_matrix.tsv - TSV file with path names and node coverages from `odgi paths` 
- sample.gafpack.gz - GAF (graph alignment format) file for a sample compressed with pigz
- bad_samples.txt - Text file with names of paths to exclude (one per line)  
- prefix - Prefix for output files

This will generate:

- prefix/combos.tsv - Cosine similarity scores for all haplotype combinations
- prefix/best_genotype.tsv - Best matching haplotype pair and score

## Installation

Cosigt requires Go 1.16+ to install.

```
git clone https://github.com/user/cosigt
cd cosigt
go build
```

## Algorithm

For each pair of haplotypes from the graph locus:

1. Sum the coverages for the pair
2. Calculate the cosine similarity between the summed coverages and the sample alignments
3. Sort pairs by cosine similarity to find best matching genotype

## License

cosigt is released under the [MIT License](LICENSE).
