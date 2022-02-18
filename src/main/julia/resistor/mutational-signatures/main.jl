using JSON
using Printf

using LightGraphs, SimpleWeightedGraphs
using ArgParse
using FASTX
using XLSX
using BioSequences

include("mut-sig-probs.jl")
include("translation-tables.jl")
using .MutSigProbs, .TranslationTables

function parse_command_line()
	s = ArgParseSettings(description="Appends codons and mutational signature probabilities to excel sheets in columns N and O")
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
		"--excel-file"
			help = "that path to an excel file on which to append probabilities"
		"--sheets"
			help = "the sheet names within the excel file to append the codon and probability columns"
			nargs = '+'
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

function toAa(s)
	AminoAcid(AA3_1[titlecase(s)])
end

function appendColumnsToExcelFile(workbook, sheets, mutation_probs_dict, fivemer_dict)
	XLSX.openxlsx(workbook, mode="rw") do xf
		for name in sheets
			sheet = xf[name]

			for row in XLSX.eachrow(sheet)
				rowNum = XLSX.row_number(row)

				if rowNum == 1
					sheet[@sprintf "N%d" rowNum] = "sig_prob"
					sheet[@sprintf "O%s" rowNum] = "codon"
					continue
				end

				resCell = XLSX.getcell(row, 1)
				wtCell = XLSX.getcell(row, 2)
				mutCell = XLSX.getcell(row, 4)

				residue = XLSX.getdata(sheet, resCell)
				wt = toAa(XLSX.getdata(sheet, wtCell))
				mut = toAa(XLSX.getdata(sheet, mutCell))

				sheet[@sprintf "N%d" rowNum] = mutation_probs_dict[(residue, wt, mut)]
				sheet[@sprintf "O%s" rowNum] = convert(String, fivemer_dict[(residue, wt, mut)])
			end
		end
	end
end

function main()
	args = parse_command_line()
	probs = parseJsonFile(args["mut-prob"])
	sequence = readSequence(args["fasta"], args["identifier"])
	init(probs)
	mutation_prob_dict, fivemer_dict = calculateMutationalProbabilities(sequence)

	if ! haskey(args, "excel-file") || ! haskey(args, "sheets")
		exit()
	end

	appendColumnsToExcelFile(args["excel-file"], args["sheets"], mutation_prob_dict, fivemer_dict)
end

main()
