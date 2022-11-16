using JSON
using LightGraphs, SimpleWeightedGraphs
using ArgParse
using FASTX
using BioSequences
using CSV
using DataFrames
using Statistics

include("src/pareto.jl")
include("src/mut-sig-probs.jl")
include("src/translation-tables.jl")
using .Pareto, .MutSigProbs, .TranslationTables

DEBUG = false

WildTypeResidueCol = "wild-type residue"
ResidueNumCol = "residue number"
MutantResidueCol = "mutant residue"
WildTypeKStarPosCol = "wild-type K* (positive)"
MutantKStarPosCol = "mutant K* (positive)"
WildTypeKStarNegCol = "wild-type K* (negative)"
MutantKStarNegCol = "mutant K* (negative)"
SigProbCol = "signature probability"
CodonCol = "codon"
HotspotCol = "hotspot count"
RankCol = "rank"

TYPES = [
	String, 
	Int, 
	String, 
	Float64, 
	Float64, 
	Float64, 
	Float64, 
	Union{Missing, Float64}, 
	Union{Missing, String}, 
	Union{Missing, Int}, 
	Union{Missing, Int}
]

function parse_command_line()
	s = ArgParseSettings(description="Outputs a CSV file with structure- and sequence-based resistance mutants.")
	@add_arg_table! s begin
		"--mut-prob"
			help = "the path to the JSON file with three-mer mutational probabilities"
			required = true
		"--fasta"
			help = "the path to the fasta file with the cDNA sequence"
			required = true
		"--identifier"
			help = "the identifier of the sequence in the fasta to extract"
			required = true
		"--csv-file"
			help = "The path to the RESISTOR CSV file"
			required = true
		"--pareto-config"
			help = "The config file for the Pareto optimization"
			required = true
		"--debug"
			help = "Saves intermediary files for inspection"
			action = :store_true
		"--c0"
			help = "The value of c0 (see Guerin et al. 2022 for more info)"
			default = 100.0
			arg_type = Float64
	end

	return parse_args(ARGS, s)
end

function parseJsonFile(f)
	open(f, "r") do reader
		JSON.parse(reader)
	end
end

function readSequence(f, identifier)
	open(FASTA.Reader, f) do reader
		for record in reader
			if FASTA.identifier(record) == identifier
				return FASTA.sequence(record)
			end
		end
	end
end

toAa(s) = AminoAcid(AA3_1[titlecase(s)])

function addCodonAndProbs(df, mutation_probs_dict, fivemer_dict)
	for ri in 1:nrow(df)
		residue = df[!, ResidueNumCol][ri]
		wt = toAa(df[!, WildTypeResidueCol][ri])
		mut = toAa(df[!,MutantResidueCol][ri])

		df[!, SigProbCol][ri] = mutation_probs_dict[(residue, wt, mut)]
		df[!, CodonCol][ri] = convert(String, fivemer_dict[(residue, wt, mut)])
	end

	df
end

function filterNaNAnd0Prob(df)
	# this also filters out NaNs
	a = df[df[!, MutantKStarPosCol] .> 0, :]

	b = a[a[!,SigProbCol] .> 0, :]
	b
end

function computeCutoff(df, c0)
	c0 * mean(10 .^ df[!, WildTypeKStarPosCol]) / mean(10 .^ df[!, WildTypeKStarNegCol])
end

function filterByCutoff(df, cutoff)
	df[10 .^ df[!, MutantKStarPosCol] ./ 10 .^ df[!, MutantKStarNegCol] .> cutoff, :]
end

function addHotspotCount(df)
	cpy = copy(df)
	gb = groupby(cpy, ResidueNumCol)

	for k in keys(gb)
		gb[k][!, HotspotCol] .= nrow(gb[k])
	end

	cpy
end

function main()
	args = parse_command_line()
	DEBUG = args["debug"]
	probs = parseJsonFile(args["mut-prob"])
	c0 = args["c0"]
	paretoConfig = args["pareto-config"]
	fastaFile = args["fasta"]
	fastaId = args["identifier"]
	csvFile = args["csv-file"]

	sequence = readSequence(fastaFile, fastaId)
	init(probs)
	mutation_prob_dict, fivemer_dict = calculateMutationalProbabilities(sequence)

	df = CSV.read(csvFile, DataFrame, types=TYPES)
	annotated = addCodonAndProbs(df, mutation_prob_dict, fivemer_dict)
	if DEBUG
		annotated |> CSV.write("sig-probs.csv")
	end

	filtered = filterNaNAnd0Prob(annotated)
	if DEBUG
		filtered |> CSV.write("no-0s.csv")
	end

	c = computeCutoff(filtered, c0)
	if DEBUG
		println(stderr, "cutoff c is $c")
	end

	resistanceMutants = filterByCutoff(filtered, c)
	if DEBUG
		resistanceMutants |> CSV.write("without-hotspot-count.csv")
	end

	withHotspotCount = addHotspotCount(resistanceMutants)
	if DEBUG
		withHotspotCount |> CSV.write("with-hotspot-count.csv")
	end

	# The maximized columns are made negative in cpy
	cpy, colnames = Pareto.preparedata(withHotspotCount, paretoConfig)
	Pareto.calcranks(cpy, colnames)

	if DEBUG
		cpy |> CSV.write("pareto-ranks-minimization.csv")
	end

	# Copy the ranks back to the last DF without negative columns
	withHotspotCount[!, RankCol] = cpy[!, RankCol]
	withHotspotCount |> CSV.write(stdout)
end

main()
