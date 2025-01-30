package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/akamensky/argparse"
	"gonum.org/v1/gonum/stat/combin"
)

// Vector represents a float64 slice for cleaner type declarations
type Vector []float64

// GetMagnitude calculates the Euclidean magnitude of two vectors with optional mask
func GetMagnitude(A, B Vector, mask []bool) float64 {
	var ALen, BLen float64
	for i, a := range A {
		if mask == nil || mask[i] {
			ALen += a * a
			BLen += B[i] * B[i]
		}
	}
	return math.Sqrt(ALen * BLen)
}

// GetDotProduct calculates the dot product of two vectors with optional mask
func GetDotProduct(A, B Vector, mask []bool) float64 {
	var dotProduct float64
	for i, a := range A {
		if mask == nil || mask[i] {
			dotProduct += a * B[i]
		}
	}
	return dotProduct
}

// GetCosineSimilarity calculates the cosine similarity between two vectors with optional mask
func GetCosineSimilarity(A, B Vector, mask []bool) float64 {
	euclMagn := GetMagnitude(A, B, mask)
	if euclMagn > 0 {
		return GetDotProduct(A, B, mask) / euclMagn
	}
	return 0
}

// ReadMask reads a file containing boolean mask values (0/1)
func ReadMask(filename string) ([]bool, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening mask file: %w", err)
	}
	defer file.Close()

	var mask []bool
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		val := strings.TrimSpace(scanner.Text())
		if val == "1" {
			mask = append(mask, true)
		} else if val == "0" {
			mask = append(mask, false)
		} else {
			return nil, fmt.Errorf("invalid mask value: %s (must be 0 or 1)", val)
		}
	}
	return mask, scanner.Err()
}

// ReadBlacklist reads a file containing paths to exclude
func ReadBlacklist(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening blacklist file: %w", err)
	}
	defer file.Close()

	var ids []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		ids = append(ids, scanner.Text())
	}
	return ids, scanner.Err()
}


// ReadGz reads gzip-compressed TSV files containing gafpack and odgi paths data
func ReadGz(filename string) ([]string, []Vector, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, fmt.Errorf("error opening gzip file: %w", err)
	}
	defer file.Close()

	gr, err := gzip.NewReader(file)
	if err != nil {
		return nil, nil, fmt.Errorf("error creating gzip reader: %w", err)
	}
	defer gr.Close()

	cr := csv.NewReader(gr)
	cr.Comma = '\t'

	// Skip header
	if _, err := cr.Read(); err != nil {
		return nil, nil, fmt.Errorf("error reading header: %w", err)
	}

	var id []string
	var cov []Vector

	for {
		rec, err := cr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, nil, fmt.Errorf("error reading record: %w", err)
		}

		f64s := make(Vector, len(rec)-1)
		for i, s := range rec[1:] {
			f64s[i], err = strconv.ParseFloat(s, 64)
			if err != nil {
				return nil, nil, fmt.Errorf("error parsing float: %w", err)
			}
		}
		cov = append(cov, f64s)
		id = append(id, rec[0])
	}

	return id, cov, nil
}

// ReadJSON reads a JSON file containing cluster information
func ReadJSON(filename string) (map[string]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening JSON file: %w", err)
	}
	defer file.Close()

	var clstr map[string]string
	if err := json.NewDecoder(file).Decode(&clstr); err != nil {
		return nil, fmt.Errorf("error decoding JSON: %w", err)
	}
	return clstr, nil
}

// WriteResults writes the analysis results to output files
func WriteResults(m *sync.Map, keys []string, clstr map[string]string, outDir, id string) error {
	if err := os.MkdirAll(outDir, os.ModePerm); err != nil {
		return fmt.Errorf("error creating output directory: %w", err)
	}

	sortedCombosFile, err := os.Create(filepath.Join(outDir, "sorted_combos.tsv"))
	if err != nil {
		return fmt.Errorf("error creating sorted_combos.tsv: %w", err)
	}
	defer sortedCombosFile.Close()

	genotypeFile, err := os.Create(filepath.Join(outDir, "cosigt_genotype.tsv"))
	if err != nil {
		return fmt.Errorf("error creating cosigt_genotype.tsv: %w", err)
	}
	defer genotypeFile.Close()

	sortedCombosWriter := csv.NewWriter(sortedCombosFile)
	sortedCombosWriter.Comma = '\t'
	defer sortedCombosWriter.Flush()

	genotypeWriter := csv.NewWriter(genotypeFile)
	genotypeWriter.Comma = '\t'
	defer genotypeWriter.Flush()

	// Write headers
	if err := genotypeWriter.Write([]string{"#sample.id", "haplotype.1", "haplotype.2", "cluster.1", "cluster.2", "cosine.similarity"}); err != nil {
		return fmt.Errorf("error writing genotype header: %w", err)
	}
	if err := sortedCombosWriter.Write([]string{"#haplotype.1", "haplotype.2", "cluster.1", "cluster.2", "cosine.similarity"}); err != nil {
		return fmt.Errorf("error writing sorted combos header: %w", err)
	}

	for i, k := range keys {
		haps := strings.Split(k, "$")
		val, _ := m.Load(k) // Fetch value from sync.Map
		cosineSimilarity := val.(float64)

		if i == 0 {
			if err := genotypeWriter.Write([]string{id, haps[0], haps[1], clstr[haps[0]], clstr[haps[1]], fmt.Sprintf("%.16f", cosineSimilarity)}); err != nil {
				return fmt.Errorf("error writing genotype data: %w", err)
			}
		}
		if err := sortedCombosWriter.Write([]string{haps[0], haps[1], clstr[haps[0]], clstr[haps[1]], fmt.Sprintf("%.16f", cosineSimilarity)}); err != nil {
			return fmt.Errorf("error writing sorted combos data: %w", err)
		}
	}

	return nil
}

// SliceContains checks if a string is contained in any string of a slice
func SliceContains(s string, ids []string) bool {
	for _, x := range ids {
		if strings.Contains(s, x) {
			return true
		}
	}
	return false
}

// SumSlices adds two float64 slices element-wise
func SumSlices(a, b Vector) Vector {
	c := make(Vector, len(a))
	for i := range a {
		c[i] = a[i] + b[i]
	}
	return c
}

func main() {
	parser := argparse.NewParser("cosigt", "genotyping loci in pangenome graphs using cosine distance")
	p := parser.String("p", "paths", &argparse.Options{Required: true, Help: "gzip-compressed tsv file with path names and node coverages from odgi paths"})
	g := parser.String("g", "gaf", &argparse.Options{Required: true, Help: "gzip-compressed gaf (graph alignment format) file for a sample from gafpack"})
	b := parser.String("b", "blacklist", &argparse.Options{Required: false, Help: "txt file with names of paths to exclude (one per line)"})
	c := parser.String("c", "cluster", &argparse.Options{Required: false, Help: "cluster json file as generated with cluster.r"})
	maskFile := parser.String("m", "mask", &argparse.Options{Required: false, Help: "boolean mask to ignore node coverages"})
	o := parser.String("o", "output", &argparse.Options{Required: true, Help: "folder prefix for output files"})
	i := parser.String("i", "id", &argparse.Options{Required: true, Help: "sample name"})

	if err := parser.Parse(os.Args); err != nil {
		log.Fatalf("Error parsing arguments: %v", err)
	}

	hapid, gcov, err := ReadGz(*p)
	if err != nil {
		log.Fatalf("Error reading paths file: %v", err)
	}

	_, bcov, err := ReadGz(*g)
	if err != nil {
		log.Fatalf("Error reading gaf file: %v", err)
	}

	var blck []string
	if *b != "" {
		blck, err = ReadBlacklist(*b)
		if err != nil {
			log.Fatalf("Error reading blacklist: %v", err)
		}
	}

	clstr, err := ReadJSON(*c)
	if err != nil {
		log.Fatalf("Error reading cluster file: %v", err)
	}

	var mask []bool
	if *maskFile != "" {
		mask, err = ReadMask(*maskFile)
		if err != nil {
			log.Fatalf("Error reading mask file: %v", err)
		}
		if len(mask) != len(gcov[0]) {
			log.Fatalf("Mask length (%d) does not match vector length (%d)", len(mask), len(gcov[0]))
		}
	}

	var seen sync.Map      // Changed to use sync.Map for concurrency safety
	results := &sync.Map{} // Changed to use sync.Map for results
	n := len(hapid)
	k := 2 // fixed ploidy

	numWorkers := runtime.NumCPU()
	jobs := make(chan []int, numWorkers)
	resultsChan := make(chan map[string]float64, numWorkers)
	var wg sync.WaitGroup

	// Start worker goroutines
	for w := 0; w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for combo := range jobs {
				localResults := make(map[string]float64)
				h1, h2 := combo[0], combo[1]
				if len(blck) == 0 || (!SliceContains(hapid[h1], blck) && !SliceContains(hapid[h2], blck)) {
					sum := SumSlices(gcov[h1], gcov[h2])
					indiv := hapid[h1] + "$" + hapid[h2]
					localResults[indiv] = GetCosineSimilarity(sum, bcov[0], mask)

					_, seenH1 := seen.LoadOrStore(h1, true)
					if !seenH1 {
						sum = SumSlices(gcov[h1], gcov[h1])
						indiv = hapid[h1] + "$" + hapid[h1]
						localResults[indiv] = GetCosineSimilarity(sum, bcov[0], mask)
					}

					_, seenH2 := seen.LoadOrStore(h2, true)
					if !seenH2 {
						sum = SumSlices(gcov[h2], gcov[h2])
						indiv = hapid[h2] + "$" + hapid[h2]
						localResults[indiv] = GetCosineSimilarity(sum, bcov[0], mask)
					}
				}
				resultsChan <- localResults
			}
		}()
	}

	// Generate combinations and send jobs
	go func() {
		gen := combin.NewCombinationGenerator(n, k)
		for gen.Next() {
			combo := gen.Combination(nil)
			jobs <- combo
		}
		close(jobs)
	}()

	// Collect results
	go func() {
		for localResults := range resultsChan {
			for k, v := range localResults {
				results.Store(k, v) // Store result in sync.Map
			}
		}
		close(resultsChan)
	}()

	wg.Wait()

	// Sort results
	var keys []string
	results.Range(func(key, value interface{}) bool {
		keys = append(keys, key.(string))
		return true
	})
	sort.SliceStable(keys, func(i, j int) bool {
		valI, _ := results.Load(keys[i])
		valJ, _ := results.Load(keys[j])
		return valI.(float64) > valJ.(float64)
	})

	// Write output
	if err := WriteResults(results, keys, clstr, *o, *i); err != nil {
		log.Fatalf("Error writing results: %v", err)
	}
}
