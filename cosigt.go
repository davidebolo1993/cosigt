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

// Small epsilon value 
const epsilon = 1e-10

// GetMagnitudeWeighted calculates the weighted Euclidean magnitude of two vectors, with optional mask and weights
func GetMagnitudeWeighted(A, B Vector, mask []bool, weights []float64) float64 {
    var ALen, BLen float64
    for i, a := range A {
        if (mask == nil || mask[i]) && (weights == nil || weights[i] > 0) {
            w := 1.0
            if weights != nil {
                w = weights[i]
            }
            ALen += w * a * a
            BLen += w * B[i] * B[i]
        }
    }
    return math.Sqrt(ALen * BLen)
}

// GetDotProductWeighted calculates the weighted dot product of two vectors, with optional mask and weights
func GetDotProductWeighted(A, B Vector, mask []bool, weights []float64) float64 {
    var dotProduct float64
    for i, a := range A {
        if (mask == nil || mask[i]) && (weights == nil || weights[i] > 0) {
            w := 1.0
            if weights != nil {
                w = weights[i]
            }
            dotProduct += w * a * B[i]
        }
    }
    return dotProduct
}

// GetCosineSimilarityWeighted calculates the weighted cosine similarity between two vectors
func GetCosineSimilarityWeighted(A, B Vector, mask []bool, weights []float64) float64 {
    euclMagn := GetMagnitudeWeighted(A, B, mask, weights)
    if euclMagn > 0 {
        return GetDotProductWeighted(A, B, mask, weights) / euclMagn
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

// ReadWeights reads a file containing float weights (one per line)
func ReadWeights(filename string) ([]float64, error) {
    file, err := os.Open(filename)
    if err != nil {
        return nil, fmt.Errorf("error opening weights file: %w", err)
    }
    defer file.Close()

    var weights []float64
    scanner := bufio.NewScanner(file)
    for scanner.Scan() {
        val := strings.TrimSpace(scanner.Text())
        if val == "" {
            continue
        }
        w, err := strconv.ParseFloat(val, 64)
        if err != nil {
            return nil, fmt.Errorf("invalid weight: %s", val)
        }
        weights = append(weights, w)
    }
    return weights, scanner.Err()
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

func WriteResults(m *sync.Map, keys []string, clstr map[string]string, outDir, id string, ploidy int) error {
    if err := os.MkdirAll(outDir, os.ModePerm); err != nil {
        return fmt.Errorf("error creating output directory: %w", err)
    }

    prefix := filepath.Base(outDir)

    sortedCombosPath := filepath.Join(outDir, fmt.Sprintf("%s.sorted_combos.tsv.gz", prefix))
    sortedCombosFile, err := os.Create(sortedCombosPath)
    if err != nil {
        return fmt.Errorf("error creating %s: %w", sortedCombosPath, err)
    }
    defer sortedCombosFile.Close()

    genotypeFile, err := os.Create(filepath.Join(outDir, fmt.Sprintf("%s.cosigt_genotype.tsv", prefix)))
    if err != nil {
        return fmt.Errorf("error creating %s.genotype.tsv: %w", prefix, err)
    }
    defer genotypeFile.Close()

    gzw := gzip.NewWriter(sortedCombosFile)
    defer gzw.Close()
    sortedCombosWriter := csv.NewWriter(gzw)
    sortedCombosWriter.Comma = '\t'
    defer sortedCombosWriter.Flush()

    genotypeWriter := csv.NewWriter(genotypeFile)
    genotypeWriter.Comma = '\t'
    defer genotypeWriter.Flush()

    hasClusterInfo := len(clstr) > 0
    genotypeHeader := []string{"#sample.id"}
    sortedCombosHeader := []string{}
    for i := 1; i <= ploidy; i++ {
        genotypeHeader = append(genotypeHeader, fmt.Sprintf("haplotype.%d", i))
        sortedCombosHeader = append(sortedCombosHeader, fmt.Sprintf("haplotype.%d", i))
    }
    if hasClusterInfo {
        for i := 1; i <= ploidy; i++ {
            genotypeHeader = append(genotypeHeader, fmt.Sprintf("cluster.%d", i))
            sortedCombosHeader = append(sortedCombosHeader, fmt.Sprintf("cluster.%d", i))
        }
    }
    genotypeHeader = append(genotypeHeader, "cosine.similarity")
    sortedCombosHeader = append(sortedCombosHeader, "cosine.similarity")

    if err := genotypeWriter.Write(genotypeHeader); err != nil {
        return fmt.Errorf("error writing genotype header: %w", err)
    }
    if err := sortedCombosWriter.Write(sortedCombosHeader); err != nil {
        return fmt.Errorf("error writing sorted combos header: %w", err)
    }

    for i, k := range keys {
        haps := strings.Split(k, "$")
        val, _ := m.Load(k)
        cosineSimilarity := val.(float64)

        genotypeRow := []string{id}
        sortedCombosRow := []string{}

        for _, hap := range haps {
            genotypeRow = append(genotypeRow, hap)
            sortedCombosRow = append(sortedCombosRow, hap)
        }
        if hasClusterInfo {
            for _, hap := range haps {
                clusterVal, exists := clstr[hap]
                if !exists {
                    clusterVal = "NA"
                }
                genotypeRow = append(genotypeRow, clusterVal)
                sortedCombosRow = append(sortedCombosRow, clusterVal)
            }
        }
        genotypeRow = append(genotypeRow, fmt.Sprintf("%.16f", cosineSimilarity))
        sortedCombosRow = append(sortedCombosRow, fmt.Sprintf("%.16f", cosineSimilarity))

        if i == 0 {
            if err := genotypeWriter.Write(genotypeRow); err != nil {
                return fmt.Errorf("error writing genotype data: %w", err)
            }
        }

        if err := sortedCombosWriter.Write(sortedCombosRow); err != nil {
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
    g := parser.String("g", "gafpack", &argparse.Options{Required: true, Help: "gzip-compressed tsv file with node coverages for a sample from gafpack"})
    b := parser.String("b", "blacklist", &argparse.Options{Required: false, Help: "txt file with names of paths to exclude (one string per line)"})
    c := parser.String("c", "cluster", &argparse.Options{Required: false, Help: "cluster json file as generated with cluster.r"})
    maskFile := parser.String("m", "mask", &argparse.Options{Required: false, Help: "boolean mask to ignore node coverages (one boolean per line)"})
    w := parser.String("w", "weights", &argparse.Options{Required: false, Help: "file with node weights (one float64 per line)"})
    ploidy := parser.Int("n", "ploidy", &argparse.Options{Required: false, Help: "ploidy level (default: 2)", Default: 2})
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

    var clstr map[string]string
    if *c != "" {
        clstr, err = ReadJSON(*c)
        if err != nil {
            log.Fatalf("Error reading cluster file: %v", err)
        }
    } else {
        clstr = make(map[string]string)
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

    var weights []float64
    if *w != "" {
        weights, err = ReadWeights(*w)
        if err != nil {
            log.Fatalf("Error reading weights file: %v", err)
        }
        if len(weights) != len(gcov[0]) {
            log.Fatalf("Weights length (%d) does not match vector length (%d)", len(weights), len(gcov[0]))
        }
    }

    var seen sync.Map
    results := &sync.Map{}
    n := len(hapid)
    k := *ploidy

    numWorkers := runtime.NumCPU()
    jobs := make(chan []int, numWorkers)
    resultsChan := make(chan map[string]float64, numWorkers)
    var wg sync.WaitGroup

    for w := 0; w < numWorkers; w++ {
        wg.Add(1)
        go func() {
            defer wg.Done()
            for combo := range jobs {
                localResults := make(map[string]float64)
                skipCombo := false
                if len(blck) > 0 {
                    for _, idx := range combo {
                        if SliceContains(hapid[idx], blck) {
                            skipCombo = true
                            break
                        }
                    }
                }
                if !skipCombo {
                    var sum Vector
                    var haplotypeIDs []string
                    for _, idx := range combo {
                        if sum == nil {
                            sum = make(Vector, len(gcov[idx]))
                            copy(sum, gcov[idx])
                        } else {
                            sum = SumSlices(sum, gcov[idx])
                        }
                        haplotypeIDs = append(haplotypeIDs, hapid[idx])
                    }
                    indiv := strings.Join(haplotypeIDs, "$")
                    localResults[indiv] = GetCosineSimilarityWeighted(sum, bcov[0], mask, weights)
                    for _, idx := range combo {
                        _, seenHap := seen.LoadOrStore(idx, true)
                        if !seenHap {
                            homoSum := make(Vector, len(gcov[idx]))
                            homoIDs := make([]string, k)
                            for j := 0; j < k; j++ {
                                if j == 0 {
                                    copy(homoSum, gcov[idx])
                                } else {
                                    homoSum = SumSlices(homoSum, gcov[idx])
                                }
                                homoIDs[j] = hapid[idx]
                            }
                            homoIndiv := strings.Join(homoIDs, "$")
                            localResults[homoIndiv] = GetCosineSimilarityWeighted(homoSum, bcov[0], mask, weights)
                        }
                    }
                }
                resultsChan <- localResults
            }
        }()
    }

    go func() {
        gen := combin.NewCombinationGenerator(n, k)
        for gen.Next() {
            combo := gen.Combination(nil)
            jobs <- combo
        }
        close(jobs)
    }()

    go func() {
        for localResults := range resultsChan {
            for k, v := range localResults {
                results.Store(k, v)
            }
        }
        close(resultsChan)
    }()

    wg.Wait()

    var keys []string
    results.Range(func(key, value interface{}) bool {
        keys = append(keys, key.(string))
        return true
    })

    sort.SliceStable(keys, func(i, j int) bool {
        valI, _ := results.Load(keys[i])
        valJ, _ := results.Load(keys[j])
        
        cosI := valI.(float64)
        cosJ := valJ.(float64)
        
        // Check if effectively equal within tolerance
        if math.Abs(cosI - cosJ) < epsilon {
            return keys[i] < keys[j]
        }
        
        return cosI > cosJ
    })

    if err := WriteResults(results, keys, clstr, *o, *i, *ploidy); err != nil {
        log.Fatalf("Error writing results: %v", err)
    }
}

