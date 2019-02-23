package main

import (
	"encoding/csv"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
)

type popDetails struct {
	family int
	price  int
}

var productCode int

const maxPrice = 10000
const univSize = 16
const popSize = 40
const kitPrice = 9000          // in cents
const familyRatio = 7.0 / 10.0 // at least 7 families in 10 products
const tournSize = 2
const minFit = 1000
const factor = 500.0 // for family diversity, the more, the better
const penalty = 0.25 // for distance of sum of prices to 9000, when below
const probMutation = 0.05
const productsFileName = "products.txt"

func getInput(filename string) map[int]popDetails {
	m := make(map[int]popDetails)
	var prod, fam, price int

	b, err := ioutil.ReadFile(filename)
	if err != nil {
		fmt.Print(err)
	}

	str := string(b) // convert content to a 'string'
	//fmt.Println(str) // print the content as a 'string'

	r := csv.NewReader(strings.NewReader(str))
	s, err := r.ReadAll()

	for i := range s {
		prod, err = strconv.Atoi(s[i][0])
		if err != nil {
			fmt.Println(err)
			os.Exit(2)
		}
		fam, err = strconv.Atoi(s[i][1])
		if err != nil {
			fmt.Println(err)
			os.Exit(2)
		}
		price, err = strconv.Atoi(s[i][2])
		if err != nil {
			fmt.Println(err)
			os.Exit(2)
		}
		m[prod] = popDetails{fam, price}
	}
	return m
}

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
func chooseCodes(univ map[int]popDetails, sz int) []int {
	allCodes := make([]int, 0)
	for key := range univ {
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
func quantifyChromossome(univ map[int]popDetails, chrom []int, ucodes []int) (int, map[int]int) {
	price := 0
	families := make(map[int]int)
	pCode := 0
	for i := 0; i < len(chrom); i++ {
		// if it's active, get price and family
		if chrom[i] == 1 {
			pCode = ucodes[i]
			u := univ[pCode]
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
	case p == 0:
		fit += 1.0
	case p < kitPrice:
		fit += 6*penalty*factor*float64(kitPrice)/float64(kitPrice-p+kitPrice) - 2*factor
	default:
		fit += 6 * factor * float64(kitPrice) / float64(p)
	}

	// 	fmt.Println("\u0394p: ", p, "  f: ", f, "  fit: ", fit)

	nrFamilies := len(f)
	for _, v := range f {
		nrProducts += v
	}

	if nrProducts != 0 {
		fit += factor * float64(nrFamilies) / float64(nrProducts)
		fit = math.Round(fit*100) / 100
	} else {
		fit += 1.0
	}
	return fit
}

func evaluatePopulation(univ map[int]popDetails, pop [][]int, ucodes []int) []float64 {
	fit := make([]float64, 0)
	p := 0
	fam := make(map[int]int)
	for i := 0; i < len(pop); i++ {
		p, fam = quantifyChromossome(univ, pop[i], ucodes)
		fit = append(fit, calculateFitness(p, fam))
		fmt.Println("Gene: ", pop[i], "   fit: ", fit[len(fit)-1], "  nrProducts: ", p, " nrFamilies: ", fam)
	}
	return fit
}

// getBestIndividual returns the individual with the best fitness
func getBestIndividual(pop [][]int, fit []float64) []int {
	maxFit := -10e9
	indxMax := -1
	best := make([]int, len(pop[0]))
	for i, v := range fit {
		if v > maxFit {
			maxFit = v
			indxMax = i
		}
	}
	copy(best, pop[indxMax])
	return best
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
		//fmt.Println("rnd1: ", rnd1, "   rnd2: ", rnd2, "   fpar1: ", fpar1, "  fpar2: ", fpar2)
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

func mutation(pop [][]int) {
	nrGenes := univSize * popSize
	nrMutations := int(probMutation * float64(nrGenes))
	fmt.Println("Doing nrMutations: ", nrMutations)
	for i := 0; i < nrMutations; i++ {
		individual := rndInt(0, popSize)
		alele := rndInt(0, univSize)
		old := pop[individual][alele]
		fmt.Println("individual: ", individual, "   alele: ", alele, "   actual: ", old)
		switch old {
		case 0:
			pop[individual][alele] = 1
		case 1:
			pop[individual][alele] = 0
		}
	}
}

func equalSlices(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}

	for i, v := range a {
		if v != b[i] {
			return false
		}
	}
	return true
}

// insertBestIndivual inserts the best chromossome of the last generation
// at the last place of the newly crossed and mutated children slice.
func insertBestIndivual(pop [][]int, best []int, fit []float64) {
	alreadyPresent := false
	var idxWorst int
	minFit := 10e9
	for i, f := range fit {
		if f < minFit {
			minFit = f
			idxWorst = i
		}
	}
	b := make([]int, len(best))
	copy(b, best)
	for _, p := range pop {
		if equalSlices(p, b) {
			alreadyPresent = true
			break
		}
	}
	if !alreadyPresent {
		pop[idxWorst] = b
	}
}

// decodePopulation
func decodePopulation(univ map[int]popDetails, pop [][]int, univCodes []int, fit []float64) {
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

		p, fam := quantifyChromossome(univ, ind, univCodes)
		fmt.Println(p, " _ ", fam)
	}
}

// averageFit
func averageFit(fit []float64) float64 {
	var total float64
	for i := 0; i < len(fit); i++ {
		total += fit[i]
	}
	return math.Round(total/float64(len(fit))*100) / 100
}

func saveProducts(m map[int]popDetails, filename string) {
	file, err := os.Create(filename)
	checkError("Cannot create file ", err)
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	for k, v := range m {
		s := []string{strconv.Itoa(k), strconv.Itoa(v.family), strconv.Itoa(v.price)}
		err := writer.Write(s)
		checkError("Cannot write to file ", err)
	}
}

func checkError(message string, err error) {
	if err != nil {
		log.Fatal(message, err)
	}
}

func main() {
	// maxPrice expressed in cents
	universe := make(map[int]popDetails)
	var elitism = true
	var bestIndividual []int

	fmt.Println()
	universe = mapUniverse(univSize)
	//saveProducts(universe, productsFileName)
	//universe = getInput(productsFileName )

	fmt.Println("Universe: ", universe)
	univCodes := chooseCodes(universe, univSize)
	fmt.Println("univCodes: ", univCodes)
	var population [][]int
	population = buildInitialPopulation(popSize, univCodes)
	fmt.Println("Population: ", population)
	fit := make([]float64, 0)
	fit_sorted := make([]float64, len(population))
	pool := make([]int, 0)
	statMax := make([]float64, 0)
	statAvg := make([]float64, 0)
	generations := 0
	for {
		fmt.Println("------------------------------------------------------------- generations: ", generations)
		fit = evaluatePopulation(universe, population, univCodes)
		if elitism {
			bestIndividual = getBestIndividual(population, fit)
		}
		copy(fit_sorted, fit)
		sort.Float64s(fit_sorted)
		fmt.Println("fit_sorted: ", fit_sorted)
		statMax = append(statMax, fit_sorted[len(fit_sorted)-1])
		statAvg = append(statAvg, averageFit(fit_sorted))
		if generations == 1 {
			break
		}

		pool = selectParents(population, fit, univCodes)
		generations += 1
		fmt.Println("pool: ", pool)
		population = crossover(population, pool)
		mutation(population)
		fmt.Println("offspring before elitism", population)
		if elitism == true {
			insertBestIndivual(population, bestIndividual, fit)
		}
		fmt.Println("offspring after elitism", population)
	}
	fmt.Println("statMax: ", statMax)
	fmt.Println("statAvg: ", statAvg)
	decodePopulation(universe, population, univCodes, fit)
}
