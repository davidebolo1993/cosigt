package main

import (
	"bufio"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
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

	/*
	Calculate magnitude
	*/

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

	/*
	Calculate dot product
	*/

	var dot_product float64

	for i := 0; i < len(A); i++ {

		dot_product += (A[i] * B[i])

	}

	return dot_product
}


func GetCosineSimilarity(A []float64, B []float64) float64 {

	/*
	Calculate cosine similarity
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
	Read file with paths to exclude - can be empty
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
	Read gafpack and odgi paths files
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

func ReadJson(s string) map[string]string {

	/*
	Read cluster json file
	*/

	clstr:=make(map[string]string)
	jsonFile,_ := os.Open(s)
	defer jsonFile.Close()
	byteValue, _ := ioutil.ReadAll(jsonFile)
	json.Unmarshal([]byte(byteValue), &clstr)

	return clstr

}




func WriteResults(m map[string]float64, keys *[]string, clstr map[string]string, s string, id string) {

	/*
	Write results to sorted_combos and cosigt_genotype
	*/

	f1, _ := os.Create(s + "/sorted_combos.tsv")
	defer f1.Close()

	f2, _ := os.Create(s + "/cosigt_genotype.tsv")
	defer f2.Close()

	w := csv.NewWriter(f1)
	w.Comma = '\t'
	defer w.Flush()

	x := csv.NewWriter(f2)
	x.Comma = '\t'
	defer x.Flush()

	for i,k:=range (*keys) {

		haps:=strings.Split(k,"$")

		if i == 0 {

			//k, fmt.Sprintf("%.16f", v)
			_ = x.Write([]string{"id", "h1", "h2", "c1", "c2", "cs"})
			_ = w.Write([]string{"h1", "h2", "c1", "c2", "cs"})
			_ = x.Write([]string{id,haps[0],haps[1],clstr[haps[0]],clstr[haps[1]],fmt.Sprintf("%.16f", m[k])})

		}

		_ = w.Write([]string{haps[0],haps[1],clstr[haps[0]],clstr[haps[1]],fmt.Sprintf("%.16f", m[k])})

	}

}

func SliceContains(s string, ids []string) bool {

	/*
	Check if substring in any string
	*/

	for _, x := range ids {

		if strings.Contains(s, x) {

			return true

		}

	}

	return false

}

func SumSlices(a []float64, b []float64) []float64 {
	/*
	Sum slices of floats
	*/

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
	c := parser.String("c", "cluster", &argparse.Options{Required: false, Help: "cluster json file as generated with cluster.py"})
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
	//read blacklist
	blck := ReadBlacklist(*b)
	//read cluster .json file
	clstr:= ReadJson(*c)
	//trace combinations we have already seen
	seen:=make(map[int]bool)
	//store results in map
	m := make(map[string]float64)
	//generate all possible diploid combination
	n := len(hapid)
	k := 2 //fixed ploidy
	gen := combin.NewCombinationGenerator(n, k)

	for gen.Next() {

		h1, h2 := gen.Combination(nil)[0], gen.Combination(nil)[1]

		if len(blck) == 0 || (!SliceContains(hapid[h1], blck) && !SliceContains(hapid[h2], blck)) { //nothing to blacklist or both ids not in blacklist

			sum := SumSlices(gcov[h1], gcov[h2])
			indiv := (hapid[h1] + "$" + hapid[h2])
			m[indiv] = GetCosineSimilarity(sum, bcov[0])

			_,ok:=seen[h1]

			if !ok {

				sum=SumSlices(gcov[h1], gcov[h1])
				indiv= (hapid[h1] + "$" + hapid[h1])
				m[indiv] = GetCosineSimilarity(sum, bcov[0])
				seen[h1] = true

			}

			_,ok=seen[h2]

			if !ok {

				sum=SumSlices(gcov[h2], gcov[h2])
				indiv= (hapid[h2] + "$" + hapid[h2])
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

	//write genotype
	outpath:=path.Clean(*o)
	_ = os.Mkdir(outpath, os.ModePerm)
	WriteResults(m,&keys,clstr,outpath,smpl[0])
}
