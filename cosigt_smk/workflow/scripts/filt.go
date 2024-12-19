package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
)

// Result represents the distance metrics between two groups
type Result struct {
	GroupA         string
	GroupB         string
	EuclideanDist  float64
	JaccardDist    float64
	CosineDissim   float64
	ManhattanDist  float64
}

func main() {
	// Configure logging to write to stderr
	log.SetOutput(os.Stderr)

	// Validate command-line arguments
	if len(os.Args) < 3 {
		log.Fatalf("Usage: %s <input_file.gz> <node_length_file>", os.Args[0])
	}

	inputFile := os.Args[1]
	nodeLengthFile := os.Args[2]

	// Read input data
	df, err := readGzippedTSV(inputFile)
	if err != nil {
		log.Fatalf("Error reading input file: %v", err)
	}

	// Read node lengths
	nodelengths, err := readNodeLengths(nodeLengthFile)
	if err != nil {
		log.Fatalf("Error reading node length file: %v", err)
	}

	// Filter nodes based on coverage differences
	filtered := filterNodes(df, nodelengths)

	// Generate group combinations
	combinations := generateCombinations(filtered["path.name"])

	// Compute distances concurrently
	results := computeDistances(filtered, combinations)

	// Track printed groups to add self-combinations
	printedGroups := make(map[string]bool)

	// Print header
	fmt.Println("group.a\tgroup.b\teuclidean.dist\tjaccard.dist\tcosine.dissim\tmanhattan.dist")

	// Print results with self-combinations
	for _, res := range results {
		// Add self-combinations with zero distances
		for _, group := range []string{res.GroupA, res.GroupB} {
			if !printedGroups[group] {
				log.Printf("Adding self-combination for group: %s", group)
				fmt.Printf("%s\t%s\t0.0000\t0.0000\t0.0000\t0.0000\n", group, group)
				printedGroups[group] = true
			}
		}

		// Print group pair results in both directions
		fmt.Printf("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n", res.GroupA, res.GroupB, res.EuclideanDist, res.JaccardDist, res.CosineDissim, res.ManhattanDist)
		fmt.Printf("%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n", res.GroupB, res.GroupA, res.EuclideanDist, res.JaccardDist, res.CosineDissim, res.ManhattanDist)
	}
}

// readGzippedTSV reads a gzipped TSV file into a map of string slices
func readGzippedTSV(filename string) (map[string][]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	gzReader, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	defer gzReader.Close()

	reader := bufio.NewReader(gzReader)

	var headers []string
	data := make(map[string][]string)
	lineNum := 0

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}

		line = strings.TrimSpace(line)
		columns := strings.Split(line, "\t")

		if lineNum == 0 {
			// Capture headers
			headers = columns
			for _, header := range headers {
				data[header] = []string{}
			}
		} else {
			// Store data
			for i, col := range columns {
				data[headers[i]] = append(data[headers[i]], col)
			}
		}
		lineNum++
	}

	return data, nil
}

// readNodeLengths reads node length file and returns a set of valid nodes
func readNodeLengths(filename string) (map[string]bool, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open node length file: %w", err)
	}
	defer file.Close()

	reader := bufio.NewReader(file)
	data := make(map[string]bool)

	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, fmt.Errorf("error reading node length line: %w", err)
		}

		columns := strings.Fields(line)
		if len(columns) < 2 {
			log.Printf("Skipping invalid line: %s", line)
			continue
		}

		nodeName := columns[0]
		length, err := strconv.Atoi(columns[1])
		if err != nil {
			log.Printf("Invalid length for node %s: %v", nodeName, err)
			continue
		}

		if length > 1 {
			data[nodeName] = true
			log.Printf("Valid node: %s (length: %d)", nodeName, length)
		}
	}

	return data, nil
}

// filterNodes filters nodes based on differences in coverage and node lengths
func filterNodes(data map[string][]string, validNodes map[string]bool) map[string][]string {
	filtered := make(map[string][]string)
	filtered["path.name"] = data["path.name"]

	for id, values := range data {
		if id == "path.name" {
			continue
		}

		// Check node length validity
		if _, ok := validNodes[id]; !ok {
			log.Printf("Discarded node %s: too small", id)
			continue
		}

		// Check for differences in the column
		var diffExists bool
		for i := 1; i < len(values); i++ {
			if values[i] != values[i-1] {
				diffExists = true
				break
			}
		}

		if diffExists {
			filtered[id] = values
		} else {
			log.Printf("Discarded node %s: no difference", id)
		}
	}

	return filtered
}

// generateCombinations creates all unique group pair combinations
func generateCombinations(groups []string) [][2]string {
	var combinations [][2]string
	for i := 0; i < len(groups); i++ {
		for j := i + 1; j < len(groups); j++ {
			combinations = append(combinations, [2]string{groups[i], groups[j]})
		}
	}
	return combinations
}

// computeDistances calculates various distance metrics concurrently
func computeDistances(data map[string][]string, combinations [][2]string) []Result {
	results := make([]Result, 0, len(combinations))
	var mu sync.Mutex
	var wg sync.WaitGroup

	for _, comb := range combinations {
		wg.Add(1)
		go func(pair [2]string) {
			defer wg.Done()

			groupA, groupB := pair[0], pair[1]
			h1 := extractValues(data, groupA)
			h2 := extractValues(data, groupB)

			// Parse values into floats
			v1 := make([]float64, len(h1))
			v2 := make([]float64, len(h2))
			for i := range h1 {
				v1[i], _ = strconv.ParseFloat(h1[i], 64)
				v2[i], _ = strconv.ParseFloat(h2[i], 64)
			}

			// Compute Euclidean distance
			eucDist := computeEuclideanDistance(v1, v2)

			// Compute Jaccard distance
			jaccDist := computeJaccardDistance(h1, h2)

			// Compute Cosine Dissimilarity
			cosineDissim := computeCosineDissimilarity(v1, v2)

			// Compute Manhattan distance
			manhattanDist := computeManhattanDistance(v1, v2)

			// Add result to the list
			res := Result{
				GroupA:         groupA,
				GroupB:         groupB,
				EuclideanDist:  eucDist,
				JaccardDist:    jaccDist,
				CosineDissim:   cosineDissim,
				ManhattanDist:  manhattanDist,
			}

			mu.Lock()
			results = append(results, res)
			mu.Unlock()
		}(comb)
	}

	wg.Wait()
	return results
}

// extractValues retrieves values for a specific group
func extractValues(data map[string][]string, group string) []string {
	// Find the row index where "path.name" matches the group
	rowIndex := -1
	for i, name := range data["path.name"] {
		if strings.TrimSpace(name) == strings.TrimSpace(group) {
			rowIndex = i
			break
		}
	}

	if rowIndex == -1 {
		log.Fatalf("Group %s not found in path.name", group)
	}

	// Extract and sort node columns
	var nodeColumns []string
	for key := range data {
		if strings.HasPrefix(key, "node.") {
			nodeColumns = append(nodeColumns, key)
		}
	}

	sort.Slice(nodeColumns, func(i, j int) bool {
		numI, _ := strconv.Atoi(strings.TrimPrefix(nodeColumns[i], "node."))
		numJ, _ := strconv.Atoi(strings.TrimPrefix(nodeColumns[j], "node."))
		return numI < numJ
	})

	// Extract values for sorted node columns
	var values []string
	for _, key := range nodeColumns {
		values = append(values, data[key][rowIndex])
	}

	return values
}

// Computation helper functions
func computeEuclideanDistance(v1, v2 []float64) float64 {
	var sumSquared float64
	for i := range v1 {
		diff := v1[i] - v2[i]
		sumSquared += diff * diff
	}
	return math.Sqrt(sumSquared)
}

func computeJaccardDistance(h1, h2 []string) float64 {
	intersection := 0
	union := len(h1) + len(h2)
	for i := range h1 {
		if h1[i] == h2[i] {
			intersection++
		}
	}
	return 1.0 - float64(intersection)/float64(union-intersection)
}

func computeCosineDissimilarity(v1, v2 []float64) float64 {
	dotProduct, magA, magB := 0.0, 0.0, 0.0
	for i := range v1 {
		dotProduct += v1[i] * v2[i]
		magA += v1[i] * v1[i]
		magB += v2[i] * v2[i]
	}
	magA, magB = math.Sqrt(magA), math.Sqrt(magB)
	cosineSim := dotProduct / (magA * magB)
	return 1 - cosineSim
}

func computeManhattanDistance(v1, v2 []float64) float64 {
	var sumAbsDiff float64
	for i := range v1 {
		sumAbsDiff += math.Abs(v1[i] - v2[i])
	}
	return sumAbsDiff
}