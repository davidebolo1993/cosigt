package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/stat/combin"
)

//all these functions assume no errors in input files

func GetCosineSimilarity(A []float64, B []float64) float64 {

	var eucl_magn, dot_product, similarity float64

	eucl_magn = GetMagnitude(A, B)

	if eucl_magn > 0 {

		dot_product = GetDotProduct(A, B)
		similarity = dot_product / eucl_magn

	} else {

		similarity = 0

	}

	return similarity
}

func GetMagnitude(A []float64, B []float64) float64 {

	var A_len, B_len float64
	A_len = 0
	B_len = 0

	for i := 0; i < len(A); i += 1 {

		A_len += (A[i] * A[i])
		B_len += (B[i] * B[i])

	}

	return math.Sqrt(A_len * B_len)
}

func GetDotProduct(A []float64, B []float64) float64 {

	var dot_product float64

	for i := 0; i < len(A); i++ {

		dot_product += (A[i] * B[i])

	}

	return dot_product
}

func ReadBlacklist(s string) []string {

	ids := make([]string, 0, 10) //don't want to re-allocate too many times

	f, _ := os.Open(s)
	defer f.Close()

	b := bufio.NewScanner(f)

	for b.Scan() {

		ids = append(ids, b.Text())

	}

	return ids

}

func ReadGz(s string) ([]string, [][]float64) {

	cov := make([][]float64, 0, 1000)
	id := make([]string, 0, 1000)

	f, _ := os.Open(s)
	defer f.Close()

	gr, _ := gzip.NewReader(f)
	defer gr.Close()

	cr := csv.NewReader(gr)
	cr.Comma = '\t'

	//skip header
	cr.Read()

	for {

		rec, err := cr.Read()

		if err == io.EOF {

			break

		}

		f64s := make([]float64, len(rec)-1) //rec[0] is the haplotype name in z.paths and sample name in x.gafpack

		for i, s := range rec[1:] {

			f64s[i], _ = strconv.ParseFloat(s, 64)

		}

		cov = append(cov, f64s)
		id = append(id, rec[0])

	}
	return id, cov
}

func WriteMap(m map[string]float64, s string) {

	f, _ := os.Create(s)
	defer f.Close()

	w := csv.NewWriter(f)
	w.Comma = '\t'
	defer w.Flush()

	for k, v := range m {

		_ = w.Write([]string{k, fmt.Sprintf("%.16f", v)})

	}

}

func SliceContains(s string, ids []string) bool {

	for _, x := range ids {

		if strings.Contains(s, x) {

			return true

		}

	}

	return false

}

func SumSlices(a []float64, b []float64) []float64 {

	c := make([]float64, len(a))

	for i := 0; i < len(a); i++ {

		c[i] = a[i] + b[i]

	}

	return c

}

func main() {

	//I/O files

	//read first table
	hapid, gcov := ReadGz(os.Args[1])
	//read second table
	smpl, bcov := ReadGz(os.Args[2])
	//read blacklist - it can be empty or not
	blck := ReadBlacklist(os.Args[3])

	//store results in map
	m := make(map[string]float64)

	n := len(hapid)
	k := 2
	gen := combin.NewCombinationGenerator(n, k)

	for gen.Next() {

		h1, h2 := gen.Combination(nil)[0], gen.Combination(nil)[1]

		if len(blck) == 0 || (!SliceContains(hapid[h1], blck) && !SliceContains(hapid[h2], blck)) { //nothing to blacklist or both ids not in blacklist

			sum := SumSlices(gcov[h1], gcov[h2])
			indiv := (hapid[h1] + "-" + hapid[h2])
			m[indiv] = GetCosineSimilarity(sum, bcov[0])

		}

	}

	//sort map
	keys := make([]string, 0, len(m))

	for key := range m {

		keys = append(keys, key)

	}

	sort.SliceStable(keys, func(i, j int) bool {

		return m[keys[i]] > m[keys[j]]

	})

	//write combinations
	_ = os.Mkdir(os.Args[4], os.ModePerm)
	WriteMap(m, os.Args[4]+"/combos.tsv")

	//write best score
	f, _ := os.Create(os.Args[4] + "/best_genotype.tsv")
	defer f.Close()
	result := fmt.Sprintf("#sample\tbest_genotype\tbest_score\n%s\t%s\t%.16f\n", smpl[0], keys[0], m[keys[0]])
	_, _ = f.WriteString(result)
}
