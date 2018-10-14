package main

import (
	"cryptopepe.io/cryptopepe-reader/pepe"
	"cryptopepe.io/cryptopepe-svg/builder/look"
	"fmt"
	"encoding/json"
	"io/ioutil"
	"bytes"
)

func main() {
	propSpec := buildTypePartDnaTable()
	buf := new(bytes.Buffer)
	enc := json.NewEncoder(buf)
	enc.SetEscapeHTML(false)
	if err := enc.Encode(propSpec); err != nil {
		panic(err)
	}
	if err := ioutil.WriteFile("output.json", buf.Bytes(), 0644); err != nil {
		panic(err)
	}
}

type Locus struct {
	Chromosome uint
	Start uint
	Len uint
}

type ExpressorAssignment struct {
	expressor pepe.GeneExpressor
	expressedLookVar func(look *look.PepeLook)*string
}

// For future simulations?
//var colorGenes = map[Locus]ExpressorAssignment{
//	Locus{0, 0, 10}:  {pepe.SkinColorExpressor, func(look *look.PepeLook) *string {return &look.Skin.Color}},
//	Locus{0, 10, 10}: {pepe.EyesColorExpressor, func(look *look.PepeLook) *string {return &look.Head.Eyes.EyeColor}},
//	Locus{0, 44, 5}:  {pepe.HeadHatColorExpressor, func(look *look.PepeLook) *string {return &look.Head.Hair.HatColor}},
//	Locus{0, 49, 5}:  {pepe.HeadHatColor2Expressor, func(look *look.PepeLook) *string {return &look.Head.Hair.HatColor2}},
//	Locus{0, 54, 6}:  {pepe.HeadHairColorExpressor, func(look *look.PepeLook) *string {return &look.Head.Hair.HatColor}},
//	Locus{1,16, 10}:  {pepe.BodyShirtColorExpressor, func(look *look.PepeLook) *string {return &look.Body.Shirt.ShirtColor}},
//	Locus{1,34, 10}:  {pepe.GlassesPrimaryColorExpressor, func(look *look.PepeLook) *string {return &look.Extra.Glasses.PrimaryColor}},
//	// yes, there are two of them, mistakes were made, now it is a feature, double-dominant alleles for secondary glasses color
//	Locus{1,44, 10}:  {pepe.GlassesSecondaryColorExpressor, func(look *look.PepeLook) *string {return &look.Extra.Glasses.SecondaryColor}},
//	Locus{1,54, 10}:  {pepe.GlassesSecondaryColorExpressor, func(look *look.PepeLook) *string {return &look.Extra.Glasses.SecondaryColor}},
//	Locus{1,64, 10}:  {pepe.BackgroundColorExpressor, func(look *look.PepeLook) *string {return &look.BackgroundColor}},
//}

var typeGenes = map[Locus]ExpressorAssignment{
	Locus{0, 20, 12}: {pepe.EyesTypeExpressor, func(look *look.PepeLook) *string {return &look.Head.Eyes.EyeType}},
	Locus{0, 32, 12}: {pepe.HeadHairTypeExpressor, func(look *look.PepeLook) *string {return &look.Head.Hair.HairType}},
	Locus{0, 60, 12}: {pepe.HeadMouthExpressor, func(look *look.PepeLook) *string {return &look.Head.Mouth}},
	Locus{1, 0, 8}:   {pepe.BodyNeckExpressor, func(look *look.PepeLook) *string {return &look.Body.Neck}},
	Locus{1, 8, 8}:   {pepe.BodyShirtTypeExpressor, func(look *look.PepeLook) *string {return &look.Body.Shirt.ShirtType}},
	Locus{1, 26, 8}:  {pepe.GlassesTypeExpressor, func(look *look.PepeLook) *string {return &look.Extra.Glasses.GlassesType}},
}

type DnaBitSpec struct {
	BitIndex uint `json:"geneIndex"`
	Chances0 map[string]float32 `json:"c0"`
	Chances1 map[string]float32 `json:"c1"`
	Chances0Dom map[string]float32 `json:"c0dom"`
	Chances1Dom map[string]float32 `json:"c1dom"`
}


// Note: DNA is read left to right (because word-like positioning, reading direction)
//    But bit indexes within genes are read right to left! (because variable length genes).
//    And chromosome (128 bits) 1 is left, chromosome 0 is right
var geneBitMasks = [12]uint32{
	1 << 0, 1 << 1, 1 << 2, 1 << 3, 1 << 4, 1 << 5, 1 << 6, 1 << 7, 1 << 8, 1 << 9, 1 << 10, 1 << 11,
}

func runTypeGeneSimulation(locus Locus, expressor pepe.GeneExpressor, lookObj *look.PepeLook, prop *string) []*DnaBitSpec {
	specs := make([]*DnaBitSpec, locus.Len, locus.Len)
	for i := uint(0); i < locus.Len; i++ {
		bitspec := new(DnaBitSpec)
		bitspec.BitIndex = i
		bitspec.Chances0 = make(map[string]float32)
		bitspec.Chances1 = make(map[string]float32)
		bitspec.Chances0Dom = make(map[string]float32)
		bitspec.Chances1Dom = make(map[string]float32)
		specs[i] = bitspec
	}
	maxAllelNum := uint32(1) << locus.Len
	for allel := uint32(0); allel < maxAllelNum; allel++ {
		// Reset property so that non-dominant values will be shown too.
		*prop = ""
		expressor(allel, lookObj)
		for i := uint(0); i < locus.Len; i++ {
			// Every value is 0 by default, just add 1
			if allel & geneBitMasks[i] == 0 {
				specs[i].Chances0[*prop] += 1.0
			} else {
				specs[i].Chances1[*prop] += 1.0
			}
		}
		// Now change the property to something already set, so that only dominant alleles will get a chance.
		// "X" will be the chance of the allel having a recessive expression.
		*prop = "X"
		expressor(allel, lookObj)
		for i := uint(0); i < locus.Len; i++ {
			// Every value is 0 by default, just add 1
			if allel & geneBitMasks[i] == 0 {
				specs[i].Chances0Dom[*prop] += 1.0
			} else {
				specs[i].Chances1Dom[*prop] += 1.0
			}
		}
	}
	//normalizer := float32(maxAllelNum)
	for _, spec := range specs {
		for k, v := range spec.Chances0 {
			spec.Chances1[k] = v
			spec.Chances0[k] = v
			spec.Chances1Dom[k] = v
			spec.Chances0Dom[k] = v
		}
	}
	return specs
}

const chromosomeSize = 128

func buildTypePartDnaTable() []*DnaBitSpec {
	table := make([]*DnaBitSpec, 256)
	tempLook := new(look.PepeLook)
	for locus, expressorAssignment := range typeGenes {
		fmt.Printf("Building bit specs for chrom: %d locus: %d, len: %d\n", locus.Chromosome, locus.Start, locus.Len)
		res := runTypeGeneSimulation(locus, expressorAssignment.expressor, tempLook, expressorAssignment.expressedLookVar(tempLook))
		chromOffset := chromosomeSize * locus.Chromosome
		begin := 256 - (chromOffset + chromosomeSize - locus.Start)
		end := 256 - (chromOffset + chromosomeSize - (locus.Start + locus.Len))
		j := locus.Len - 1
		for i := begin; i < end; i++ {
			table[i] = res[j]
			j--
		}
	}
	return table
}
