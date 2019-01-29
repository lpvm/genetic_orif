package main

import (
	"fmt"
	"math/rand"
)

type popDetails struct {
	family int
	price  int
}

var productCode int
var universe = make(map[int]popDetails)

const maxPrice = 10000
const univSize = 100
const popSize = 8
const indivSize = 10 // 10 products, means 20 aleles

func rndInt(lo, hi int) int {
	return rand.Intn(hi-lo) + lo
}

func mapUniverse(nr int) map[int]popDetails {
	p := make(map[int]popDetails)
	for i := 0; i < nr; i++ {
		p[rndInt(10000, 99999)] = popDetails{family: rndInt(0, 10), price: rndInt(0, maxPrice)}
	}
	return p
}

func randEnabled() int {
	if r := rand.Float64(); r < 0.5 {
		return 0
	}
	return 1
}

// buildGene builds each gene according to this structure
// | 1 or 0 | code number | 1 or 0 | code number| ...
// | 1 or 0 | means product is enable (1) or disabled (0)
// | code number | the code for each product as received as parameter
// Each gene has 10 possible products, i.e., ten pairs of | 1 or 0| code number |
func buildGene(codes []int) []int {
	var g []int
	var rndC int
	// i := 0
	for {
		rndC = codes[rndInt(0, len(codes))]
		if !contains(g, rndC) {
			g = append(g, randEnabled())
			g = append(g, rndC)
	// 		i++
		}
		if len(g) == indivSize*2 {
			break
		}
	}
	return g
}

func contains(sl []int, n int) bool {
	for _, a := range sl {
		if a == n {
			return true
		}
	}
	return false
}

// uniqCodes builds a []int with popSize product codes that are unique
func uniqCodes(all []int) []int {
	uniq := make([]int, popSize)
	for c := 0; c < popSize; c++ {
		r := rndInt(0, univSize)
		if !contains(uniq, r) {
			uniq = append(uniq, r)
		}
		if len(uniq) == popSize {
			return uniq
		}
	}
	return uniq
}

func buildInitialPopulation(sz int, all []int) [][]int {
	var chrom [][]int
	for i := 0; i < sz; i++ {
		chrom = append(chrom, buildGene(all))
	}
	return chrom
}

// chooseCodes converts all codes that are keys of the universe
// and make a slice of them
func chooseCodes(sz int) []int {
	allCodes := make([]int, univSize)
	for key := range universe {
		allCodes = append(allCodes, key)
	}
	return allCodes
}

func main() {
	// maxPrice expressed in cents

	universe = mapUniverse(univSize)
	fmt.Println("Universe: ", universe)
	univCodes := chooseCodes(popSize)
	population := buildInitialPopulation(popSize, univCodes)
	fmt.Println("Population: ", population)

}
