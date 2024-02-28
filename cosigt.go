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
	"path"

	"github.com/akamensky/argparse"
	"gonum.org/v1/gonum/stat/combin"
)



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


func GetCosineSimilarity(A []float64, B []float64) float64 {

	/*
	Simple calculation of cosine similarity
	*/


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



func ReadBlacklist(s string) []string {

	/*
	Samples we may want to exclude, one-per-line
	*/


	ids := make([]string, 0, 100)

	f, _ := os.Open(s)
	defer f.Close()

	b := bufio.NewScanner(f)

	for b.Scan() {

		ids = append(ids, b.Text())

	}

	return ids

}

func ReadGz(s string) ([]string, [][]float64) {

	/*
	Read the gzip-compressed files
	*/

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

	/*
	Parsing of arguments with argparse
	*/

	parser := argparse.NewParser("cosigt", "genotyping loci in pangenome graphs using cosine distance")

	p := parser.String("p", "paths", &argparse.Options{Required: true, Help: "gzip-compressed tsv file with path names and node coverages from odgi paths"})
	g := parser.String("g", "gaf", &argparse.Options{Required: true, Help: "gzip-compressed gaf (graph alignment format) file for a sample from gafpack"})
	b := parser.String("b", "blacklist", &argparse.Options{Required: false, Help: "txt file with names of paths to exclude (one per line)"})
	o := parser.String("o", "output", &argparse.Options{Required: true, Help: "folder prefix for output files"})

	err := parser.Parse(os.Args)

	if err != nil  {

		fmt.Print(parser.Usage(err))
		os.Exit(1)

	}

	//read first table
	hapid, gcov := ReadGz(*p)
	//read second table
	smpl, bcov := ReadGz(*g)
	//read blacklist - it can be empty or not
	blck := ReadBlacklist(*b)
	//trace combinations we have already seen
	seen:=make(map[int]bool)
	//store results in map
	m := make(map[string]float64)
	//generate all possible diploid combination
	n := len(hapid)
	k := 2 //fixed ploidy - this can be however adjusted
	gen := combin.NewCombinationGenerator(n, k)

	for gen.Next() {

		h1, h2 := gen.Combination(nil)[0], gen.Combination(nil)[1]

		if len(blck) == 0 || (!SliceContains(hapid[h1], blck) && !SliceContains(hapid[h2], blck)) { //nothing to blacklist or both ids not in blacklist

			sum := SumSlices(gcov[h1], gcov[h2])
			indiv := (hapid[h1] + "-" + hapid[h2])
			m[indiv] = GetCosineSimilarity(sum, bcov[0])

			_,ok:=seen[h1]

			if !ok {

				sum=SumSlices(gcov[h1], gcov[h1])
				indiv= (hapid[h1] + "-" + hapid[h1])
				m[indiv] = GetCosineSimilarity(sum, bcov[0])
				seen[h1] = true

			}

			_,ok=seen[h2]

			if !ok {

				sum=SumSlices(gcov[h2], gcov[h2])
				indiv= (hapid[h2] + "-" + hapid[h2])
				m[indiv] = GetCosineSimilarity(sum, bcov[0])
				seen[h2] = true

			}


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
	outpath:=path.Clean(*o)
	_ = os.Mkdir(outpath, os.ModePerm)

	combos:=path.Clean(outpath + "/combos.tsv")
	WriteMap(m, combos)

	//write best score
	best:=path.Clean(outpath + "/best_genotype.tsv")
	f, _ := os.Create(best)
	defer f.Close()
	result := fmt.Sprintf("#sample\tbest_genotype\tbest_score\n%s\t%s\t%.16f\n", smpl[0], keys[0], m[keys[0]])
	_, _ = f.WriteString(result)	

	//all done

}
