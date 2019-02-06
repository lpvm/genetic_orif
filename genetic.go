package main

import (
	"fmt"
	"math/rand"
	"sort"
)

type popDetails struct {
	family int
	price  int
}

var productCode int
var universe = make(map[int]popDetails)

const maxPrice = 10000
const univSize = 10
const popSize = 8
const kitPrice = 9000          // in cents
const familyRatio = 7.0 / 10.0 // at least 7 families in 10 products
const tournSize = 2
const minFit = 1000

func rndInt(lo, hi int) int {
	return rand.Intn(hi-lo) + lo
}

// mapUniverse simulates the entire set of products available, with
// product code, family of product and price
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
	for c := 0; c < univSize; c++ {
		g = append(g, randEnabled())
		g = append(g, codes[c])
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
// and return sum of price of products and number of families
func quantifyGene(ch []int) (int, map[int]int) {
	active := false
	price := 0
	families := make(map[int]int)
	for i := 0; i < len(ch); i++ {
		// if it's active, set flag
		if i%2 == 0 && ch[i] == 1 {
			active = true
		}
		// if it's active, get price and family
		if i%2 == 1 && active {
			u := universe[ch[i]]
			price += u.price
			family := u.family
			families[family] += 1
			active = false
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
		fit -= penalty * float64(kitPrice-p)
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

func evaluatePopulation(pop [][]int) []float64 {
	fit := make([]float64, 0)
	p := 0
	fam := make(map[int]int)
	for i := 0; i < len(pop); i++ {
		p, fam = quantifyGene(pop[i])
		fmt.Println("Gene: ", pop[i])
		fit = append(fit, calculateFitness(p, fam))
	}
	return fit
}

func rescale(f []float64) []float64 {
	var lowest float64 = 1e300
	newFit := make([]float64, 0)
	for n := 0; n < len(f); n++ {
		if f[n] < lowest {
			lowest = f[n]
		}
	}
	addFit := minFit - lowest
	fmt.Println("adding: ", addFit, " because of lowest: ", lowest)
	for n := 0; n < len(f); n++ {
		newFit = append(newFit, f[n]+addFit)
	}
	return newFit
}

// selectParents use tournament to select parents of offspring based on their
// fitness
// returns a slice of the code numbers of parents, so that
// [par1, par1, par2, par2, par3, par3, ...]
func selectParents(p [][]int, fit []float64) []int {
	var p1 int
	var p2 int
	var fp1 float64
	var fp2 float64
	parents := make([]int, 0)
	for n := 0; n < 2*popSize; n++ {
		p1 = rndInt(0, popSize)
		p2 = rndInt(0, popSize)
		fp1 = fit[p1]
		fp2 = fit[p2]
		switch {
		case fp1 > fp2:
			parents = append(parents, p1)
		case fp2 > fp1:
			parents = append(parents, p2)
		case fp1 == fp2:
			parents = append(parents, p1)
		}
	}
	return parents
}

func crossover(pop [][]int, parents []int) [][]int {
	children := make([][]int, len(pop))
	child := make([]int, 0)
	last := -1
	indivSize := len(pop[0]) // assuming all individuals have equal length
	parent1 := make([]int, 0)
	parent2 := make([]int, 0)
	var crosspoint int
	for n := 0; n < len(parents); n++ {
		if n%2 == 0 {
			last = parents[n]
		} else {
			crosspoint = rndInt(1, indivSize)
			parent1 = pop[last]
			parent2 = pop[parents[n]]
			child = parent1[0:crosspoint]
			child = append(child, parent2[crosspoint:]...)
		}
		children[n / 2] = child
	}
	return children
}

func main() {
	// maxPrice expressed in cents

	universe = mapUniverse(univSize)
	fmt.Println("Universe: ", universe)
	univCodes := chooseCodes(univSize)
	fmt.Println("univCodes: ", univCodes)
	population := buildInitialPopulation(popSize, univCodes)
	fmt.Println("Population: ", population)
	fit := make([]float64, 0)
	pool := make([]int, 0)
	stat := make([][]float64, 0)
	statGen := make([]float64, 2)
	generations := 0
	for {
		fit = evaluatePopulation(population)
		fmt.Println("-------------------------------------------------------------")
		fmt.Println("fit1: ", fit)
		fit = rescale(fit)
		fmt.Println("fit2: ", fit)
		pool = selectParents(population, fit)
		sort.Float64s(fit)
		statGen[0] = fit[0]
		statGen[1] = fit[len(fit) - 1]
		stat = append(stat, statGen)
		fmt.Println("pool: ", pool)
		population = crossover(population, pool)
		fmt.Println("offspring ", population)
		if generations == 3 {
				break
		}
		generations += 1
	}
	fmt.Println("stat: ", stat)
}
