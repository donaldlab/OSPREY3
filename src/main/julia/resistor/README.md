Resistor
==

A program that identifies resistance mutations and ranks them using Pareto optimization.

The [Resistor algorithm](https://doi.org/10.1016/j.cels.2022.09.003) uses Pareto optimization to predict resistance mutations. 
Resistor optimizes over positive and negative design, mutation counts per position, and the probability that one amino acid will mutate into another amino acid. 
These cancer type-specific probabilities are extrapolated from Alexandrov et al's "Signatures of mutational processes in human cancer" (Nature 2013, doi: 10.1038/nature12477), in a manner described in Kaserer & Blagg's _Combining Mutational Signatures, Clonal Fitness, and Drug Affinity to Define Drug-Specific Resistane Mutations in Cancer_ (Cell Chemical Biology 2018, doi: 10.1016/j.chembiol.2018.07.013).

Once we determine the probability that a particular DNA codon will mutate to another DNA codon, we need to convert the DNA codon probabilities to amino acid probabilities.  
The way to do this with this program is to run:

```
julia --project=. main.jl \
	--mut-prob MUT-PROB \
	--fasta FASTA \
	--identifier IDENT \
	--csv-file CSV \
	--pareto-config CONFIG
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

- FASTA is a fasta formatted file with the cDNA sequence for the program you're investigating.

- IDENT is the identifer for the record in the fasta file for the cDNA sequence (since fasta files can contain more than one sequence).

- CSV is a comma separated value file in the RESISTOR table format with the positive and negative K\* scores filled out. 
A template file is available in the program distribution.

- CONFIG is a settings file in JSON format that contains the names of the columns and whether a column should be maximized or minimized. 
The content of the settings file we used for RESISTOR, for example, looked like:

{
	"columns": [
		{
			"name": "Sig Prob",
			"optimize": "maximize"
		},
		{
			"name": "K* ATP",
			"optimize": "maximize"
		},
		{
			"name": "K* Gefitinib",
			"optimize": "minimize"
		},
		{
			"name": "Count",
			"optimize": "maximize"
		}
	]
}

This indicates that there are four columns to optimize over, and that the
numerical values in the column "Sign Prob" should be maximized, "K\* ATP"
maximized, "K\* Gefitinib", and "Count" maximized.

The program then puts the mutational signature probability in column N and the codon in column O.
