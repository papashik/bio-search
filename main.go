package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

func main() {
	file, err := os.Open("uniprot_sprot.fasta")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	var sequences []string
	var currentSequence strings.Builder
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.TrimSpace(line) == "" {
			continue
		}

		if strings.HasPrefix(line, ">") {
			if currentSequence.Len() > 0 {
				sequences = append(sequences, currentSequence.String())
				currentSequence.Reset()
			}
			continue
		}
		cleanLine := strings.ReplaceAll(strings.TrimSpace(line), " ", "")
		currentSequence.WriteString(cleanLine)
	}

	if currentSequence.Len() > 0 {
		sequences = append(sequences, currentSequence.String())
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	fmt.Printf("Найдено %d последовательностей:\n", len(sequences))
}
