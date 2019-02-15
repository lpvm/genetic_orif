package main

import (
	"fmt"
	"math"
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
const univSize = 16
const popSize = 8
const kitPrice = 9000          // in cents
const familyRatio = 7.0 / 10.0 // at least 7 families in 10 products
const tournSize = 2
const minFit = 1000
const factor = 500.0 // for family diversity, the more, the better
const penalty = 2.0  // for distance of sum of prices to 9000, when below

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

// chooseCodes converts all codes that are keys of the universe
// and make a slice of them
func chooseCodes(sz int) []int {
	allCodes := make([]int, 0)
	for key := range universe {
		allCodes = append(allCodes, key)
	}
	sort.Ints(allCodes)
	return allCodes
}

// buildGene builds each gene according to this structure
// | 1 or 0 | 1 or 0 | ...
// | 1 or 0 | means product is enable (1) or disabled (0)
// each position refers to the same position in univCodes
// Each gene has 10 possible products, i.e., ten | 1 or 0|
func buildGene(codes []int) []int {
	var g []int
	for c := 0; c < univSize; c++ {
		g = append(g, randEnabled())
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

// quantifyChromossome takes a chromossome like
// [1 0 0 1 0 1 1 1 1 0]
// and return sum of price of products and number of families
func quantifyChromossome(chrom []int, ucodes []int) (int, map[int]int) {
	price := 0
	families := make(map[int]int)
	pCode := 0
	for i := 0; i < len(chrom); i++ {
		// if it's active, get price and family
		if chrom[i] == 1 {
			pCode = ucodes[i]
			u := universe[pCode]
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
	switch {
	case p < kitPrice:
		fit += 6 * penalty * factor * float64(kitPrice) / float64(kitPrice-p+kitPrice)
	default:
		fit += 6 * factor * float64(kitPrice) / float64(p)
	}

	// 	fmt.Println("\u0394p: ", p, "  f: ", f, "  fit: ", fit)

	nrFamilies := len(f)
	for _, v := range f {
		nrProducts += v
	}

	fit += factor * float64(nrFamilies) / float64(nrProducts)
	fit = math.Round(fit*100) / 100

	fmt.Println("\u0394p: ", p, "  f: ", f, "  fit: ", fit, "  nrProducts: ", nrProducts, " nrFamilies: ", nrFamilies)
	return fit
}

func evaluatePopulation(pop [][]int, ucodes []int) []float64 {
	fit := make([]float64, 0)
	p := 0
	fam := make(map[int]int)
	for i := 0; i < len(pop); i++ {
		p, fam = quantifyChromossome(pop[i], ucodes)
		fit = append(fit, calculateFitness(p, fam))
		fmt.Println("Gene: ", pop[i], "   fit: ", fit[len(fit)-1])
	}
	return fit
}

// selectParents use tournament to select parents of offspring based on their
// fitness
// returns a slice of the code numbers of parents, so that
// [par1, par1, par2, par2, par3, par3, ...]
func selectParents(p [][]int, fit []float64, ucodes []int) []int {
	var rnd1 int
	var rnd2 int
	var fpar1 float64
	var fpar2 float64
	parents := make([]int, 0)
	for n := 0; n < 2*popSize; n++ {
		rnd1 = rndInt(0, popSize)
		for {
			rnd2 = rndInt(0, popSize)
			if rnd2 != rnd1 {
				break
			}
		}
		fpar1 = fit[rnd1]
		fpar2 = fit[rnd2]
		switch {
		case fpar1 > fpar2:
			parents = append(parents, rnd1)
		case fpar2 > fpar1:
			parents = append(parents, rnd2)
		case fpar1 == fpar2:
			parents = append(parents, rnd1)
		}
		fmt.Println("rnd1: ", rnd1, "   rnd2: ", rnd2, "   fpar1: ", fpar1, "  fpar2: ", fpar2)
	}
	return parents
}

func crossover(pop [][]int, parents []int) [][]int {
	children := make([][]int, len(pop))
	child := make([]int, len(pop[0]))
	parent1 := make([]int, 0)
	parent2 := make([]int, 0)
	indivSize := len(pop[0]) // assuming all individuals have equal length
	var crosspoint int
	fmt.Println("pop: ")
	for n := 0; n < len(pop); n++ {
		fmt.Println("n: ", n, "   ", pop[n])
	}
	for n := 0; n < len(parents); n++ {
		if n%2 == 1 {
			crosspoint = rndInt(1, indivSize)
			parent1 = pop[parents[n-1]]
			parent2 = pop[parents[n]]
			child = nil
			child = make([]int, len(parent1[0:crosspoint]))
			copy(child, parent1[0:crosspoint])
			child = append(child, parent2[crosspoint:]...)
			fmt.Println("crosspoint: ", crosspoint, "   par1: ", parents[n-1], "  ", parent1, "   par2: ", parents[n], "   ", parent2, "    child: ", child)
		}
		children[n/2] = child
	}
	return children
}

func decodePopulation(pop [][]int, univCodes []int, fit []float64) {
    fmt.Println(univCodes)
	for i := 0; i < len(pop); i++ {
		ind := pop[i]
		fmt.Print(ind, " -> ")
		for p := 0; p < len(ind); p++ {
			if ind[p] == 1 {
				fmt.Print(univCodes[p], "  ")
			}

		}
		fmt.Print("\t\t", fit[i], " : ")

		p, fam := quantifyChromossome(ind, univCodes)
		fmt.Println(p, " _ ", fam)
	}
}

func main() {
	// maxPrice expressed in cents

	fmt.Println()
	universe = mapUniverse(univSize)
	fmt.Println("Universe: ", universe)
	univCodes := chooseCodes(univSize)
	fmt.Println("univCodes: ", univCodes)
	population := buildInitialPopulation(popSize, univCodes)
	fmt.Println("Population: ", population)
	fit := make([]float64, 0)
	fit_sorted := make([]float64, len(population))
	pool := make([]int, 0)
	stat := make([]float64, 0)
	generations := 0
	for {
		fmt.Println("------------------------------------------------------------- generations: ", generations)
		fit = evaluatePopulation(population, univCodes)
		fmt.Println("fit2: ", fit)
		pool = selectParents(population, fit, univCodes)
		copy(fit_sorted, fit)
		sort.Float64s(fit_sorted)
		fmt.Println("fit_sorted: ", fit_sorted)
		stat = append(stat, fit_sorted[len(fit_sorted)-1])
		fmt.Println("pool: ", pool)
		population = crossover(population, pool)
		fmt.Println("offspring ", population)
		if generations == 5 {
			break
		}
		generations += 1
	}
	fmt.Println("stat: ", stat)
	decodePopulation(population, univCodes, fit)
}
