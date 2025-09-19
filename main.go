package main

import (
	"bufio"
	"cmp"
	"fmt"
	"log"
	"math"
	"os"
	"slices"
)

const (
	TopCount = 50

	Coef1 = 1
	Coef2 = 1
)

type Seq struct {
	Name      string
	Sequence  string
	Positions map[byte][]int
	Letters   [26]int
	Result1   float64
	Result2   float64
}

func main() {
	filePath := "/Users/p00p/Downloads/uniprot_sprot.fasta"

	sequences, err := parseFastaFile(filePath)
	if err != nil {
		fmt.Printf("Error parsing file: %v\n", err)
		return
	}

	fmt.Printf("Parsed %d sequences\n", len(sequences))

	seq := sequences[0]
	fmt.Printf("Name: %s\n", seq.Name)
	fmt.Printf("Sequence length: %d\n", len(seq.Sequence))
	fmt.Printf("Letters count: %v\n", seq.Letters)
	if positions, exists := seq.Positions['A']; exists && len(positions) > 0 {
		fmt.Printf("Positions for 'A': %v\n", positions)
	}

	processPositions(seq, sequences)
	countCosinuses(seq, sequences)

	slices.SortFunc(sequences, func(a, b Seq) int {
		return cmp.Compare(
			cmp.Compare(a.Result1, b.Result1)*Coef1,
			cmp.Compare(a.Result2, b.Result2)*Coef2,
		)
	})

	for i := 0; i < TopCount; i++ {
		fmt.Printf(">%s\n", sequences[i].Name)
		fmt.Println(sequences[i].Sequence)
	}
}

func parseFastaFile(filePath string) ([]Seq, error) {
	file, err := os.Open(filePath)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := bufio.NewReader(file)
	var sequences []Seq
	var currentSeq *Seq
	var inHeader bool
	var currentHeader []byte
	var currentSequence []byte
	position := 0

	for {
		b, err := reader.ReadByte()
		if err != nil {
			if currentSeq != nil {
				currentSeq.Sequence = string(currentSequence)
				sequences = append(sequences, *currentSeq)
			}
			log.Println(err)
			break
		}

		switch {
		case b == '>':
			if currentSeq != nil {
				currentSeq.Sequence = string(currentSequence)
				sequences = append(sequences, *currentSeq)
			}

			currentSeq = &Seq{
				Positions: make(map[byte][]int),
			}
			currentHeader = make([]byte, 0, 150)
			currentSequence = make([]byte, 0, 500)
			inHeader = true
			position = 0

		case b == '\n' || b == '\r':
			if inHeader {
				if len(currentHeader) > 0 {
					currentSeq.Name = string(currentHeader)
				}
				inHeader = false
			} else if currentSeq != nil {
				continue
			}

		default:
			if inHeader {
				currentHeader = append(currentHeader, b)
			} else if currentSeq != nil {
				currentSequence = append(currentSequence, b)
				currentSeq.Letters[b-'A']++
				currentSeq.Positions[b] = append(currentSeq.Positions[b], position)
				position++
			}
		}
	}

	if currentSeq != nil {
		currentSeq.Sequence = string(currentSequence)
		sequences = append(sequences, *currentSeq)
	}

	return sequences, nil
}

func countCosinuses(seq Seq, sequences []Seq) {
	for i := range sequences {
		sequences[i].Result2 = cos(seq.Letters, sequences[i].Letters)
	}
}

func cos(v1, v2 [26]int) float64 {
	dot := 0
	mag1, mag2 := 0.0, 0.0

	for i := 0; i < 26; i++ {
		dot += v1[i] * v2[i]
		mag1 += float64(v1[i] * v1[i])
		mag2 += float64(v2[i] * v2[i])
	}

	if mag1 == 0 || mag2 == 0 {
		return 0
	}
	return float64(dot) / (math.Sqrt(mag1) * math.Sqrt(mag2))
}

func processPositions(seqInput Seq, sequences []Seq) {
	for i := range sequences {
		positions := sequences[i].Positions
		positionsInput := seqInput.Positions
		for letter, letterPos := range positions {
			letterPosInput := positionsInput[letter]
			if len(letterPosInput) > 0 {
				if len(letterPosInput) != len(letterPos) {
					pos1, pos2 := makeEqualLen(letterPosInput, letterPos)
					var sumD float64
					for i := 0; i < len(pos1); i++ {
						d := (pos1[i] - pos2[i])
						sumD += float64(d * d)
					}
					sequences[i].Result1 = sumD / float64(len(pos1))
				}
			}
		}
	}
}

func makeEqualLen(shorter, longer []int) ([]int, []int) {
	resultShort := make([]int, len(shorter))
	copy(resultShort, shorter)
	resultLong := make([]int, len(longer))
	copy(resultLong, longer)

	for len(resultShort) < len(resultLong) {
		maxDiff := 0
		insertIndex := 0

		for i := 0; i < len(resultShort)-1; i++ {
			diff := int(math.Abs(float64(resultShort[i] - resultShort[i+1])))
			if diff > maxDiff {
				maxDiff = diff
				insertIndex = i
			}
		}
		average := (resultShort[insertIndex] + resultShort[insertIndex+1]) / 2
		resultShort = insertAt(resultShort, insertIndex+1, average)
	}

	return resultShort, resultLong
}

func insertAt(arr []int, index int, value int) []int {
	if index < 0 || index > len(arr) {
		return arr
	}
	newArr := make([]int, len(arr)+1)
	copy(newArr[:index], arr[:index])
	newArr[index] = value
	copy(newArr[index+1:], arr[index:])
	return newArr
}
