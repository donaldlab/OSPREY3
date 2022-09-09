using ArgParse
using CSV
using JSON
using DataFrames

include("src/pareto.jl")
using .Pareto

csv_file = "csvfile"
config_file = "config"

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "csvfile"
            help = "The CSV file to run Pareto optimization on"
            required = true
        "config"
            help = "The config file for the Pareto optimization"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    csv = parsed_args["csvfile"]
    config = parsed_args["config"]
    df = Pareto.preparedata(csv, config)
    ranks = Pareto.calcranks(df)
    newcsv = DataFrames.hcat(CSV.read(csv, DataFrame), ranks)
    CSV.write(stdout, newcsv)
end

main()