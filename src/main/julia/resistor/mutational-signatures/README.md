Resistor - Mutational Signatures
==

A program to calculate cancer-specific mutational signature probabilities.

The [Resistor algorithm](http://textpla.in) uses Pareto optimization to predict resistance mutations. Resistor optimizes over positive and negative design, mutation counts per position, and the probability that one amino acid will mutate into another amino acid. These cancer type-specific probabilities are extrapolated from Alexandrov et al's "Signatures of mutational processes in human cancer" (Nature 2013, doi: 10.1038/nature12477), in a manner described in Kaserer & Blagg's _Combining Mutational Signatures, Clonal Fitness, and Drug Affinity to Define Drug-Specific Resistane Mutations in Cancer_ (Cell Chemical Biology 2018, doi: 10.1016/j.chembiol.2018.07.013).

Once we determine the probability that a particular DNA codon will mutate to another DNA codon, we need to convert the DNA codon probabilities to amino acid probabilities.  The way to do this with this program is to run:

```
julia --project=. main.jl \
	--mut-prob MUT-PROB \
	--fasta FASTA \
	--identifier IDENT \
	--excel-file EXCEL-FILE \
	--sheets SHEET [SHEET...]
```

- MUT-PROB is a JSON file containing the codon-to-codon probabilities, e.g:
```
{
    "ACA": {
        "AAA": 0.000889121,
        "AGA": 0.000291336,
        "ATA": 0.001601212

    },
    "ACC": {
        "AAC": 0.000904956,
        "AGC": 0.000244258,
        "ATC": 0.016219572
    },
	...
```

- EXCEL-FILE is an Excel document.

- SHEET is the name of a sheet in EXCEL-FILE. Each sheet must have the following columns:
    * In column A, the residue number.
    * In column B, the three-letter amino acid abbreviation of the wildtype
    * In column D, the three-letter amino acid abbreviation of the mutant

The program then puts the mutational signature probability in column N and the codon in column O.

- FASTA is a fasta file with the protein's coding DNA

- IDENT is the identifier for the record within the fasta file containing the cDNA.
