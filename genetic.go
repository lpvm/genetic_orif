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
const indivSize = 10 				// 10 products, means 20 aleles
const kitPrice = 9000				// in cents
const familyRatio = 7.0/10.0		// at least 7 families in 10 products

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
	for {
		rndC = codes[rndInt(0, len(codes))]
		if !contains(g, rndC) {
			g = append(g, randEnabled())
			g = append(g, rndC)
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
	allCodes := make([]int, 0)
	for key := range universe {
		allCodes = append(allCodes, key)
	}
	return allCodes
}

// quantifyGene takes a chromossome like
// [0 0 0 70026 0 31229 0 67967 1 86397 0 6327 4 1 87620 1 14976 0 98039 1 90378]
// and calculates its fitness
func quantifyGene(ch []int) (int, map[int]int) {
	active := false
	price := 0
	families := make(map[int]int)
	for i := 0; i < len(ch); i++ {
		if i%2 == 0 && ch[i] == 1 {
			active = true
		}
		if i%2 == 1 && active {
			u := universe[ch[i]]
			price += u.price
			family := u.family
			families[family] += 1
		}
	}
	return price, families
}

// calculateFitness takes the sum of prices of the products on each gene,
// in cents, and the map of nr. of products per family and calculates the
// fitness according to:
// 1. the distance from p to 9000.  Should not be lower
// 2. the product should be widespread per family
func calculateFitness(p int, f map[int]int) float64 {
	nrProducts := 0
	fit := 0.0
	factor := 500.0
	penalty := 2.0
    switch {
		case p > kitPrice:
				fit -= float64(p - kitPrice)
		case p < kitPrice:
				fit -= penalty * float64(kitPrice - p)
	}
	fmt.Println("p: ", p, "  f: ", f, "  fit: ", fit)

	nrFamilies := len(f)
	for _, v := range f {
		nrProducts += v
	}

	fit += factor * float64(nrProducts) / float64(nrFamilies)
	

	fmt.Println("p: ", p, "  f: ", f, "  fit: ", fit, "  nrProducts: ", nrProducts, " nrFamilies: ", nrFamilies)
	return fit
}

func evaluateChromossome(pop [][]int) map[int]float64 {
	fit := make(map[int]float64)
	p := 0
	f := make(map[int]int)
	for i := 0; i < len(pop); i++ {
		p, f = quantifyGene(pop[i])
		fmt.Println("Gene: ", pop[i])
		fit[i] = calculateFitness(p, f)

	}
	return fit
}

func main() {
	// maxPrice expressed in cents

	universe = mapUniverse(univSize)
	fmt.Println("Universe: ", universe)
	univCodes := chooseCodes(popSize)
	fmt.Println("univCodes: ", univCodes)
	population := buildInitialPopulation(popSize, univCodes)
	fmt.Println("Population: ", population)
	evaluateChromossome(population)

}
