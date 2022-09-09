module Pareto

using DataFrames
using JSON
using CSV

export preparedata, calcranks

RankColumn = "Rank"

"Reads the columns from the csv specified in the settings and turns all maximizes into minimizes"
function preparedata(csv, config)

    settings = JSON.parsefile(config)
    cols = settings["columns"]
    colnames = map(spec -> spec["name"], cols)
    df = CSV.read(csv, DataFrame, select=colnames)
    
    for spec in cols
        (name, optimization) = spec["name"], spec["optimize"]
        if optimization == "maximize"
            df[!, name] = -df[:, name]
        end
    end

    return df
end

"Returns an (n, 1) dimensional DataFrame corresponding to the rank of each row in `df`"
function calcranks(df)
    numrow = DataFrames.nrow(df)
    numcol = DataFrames.ncol(df)
    endrank = zeros(Int, numrow)
    dominates = [Int[] for _ in 1:numrow]
    dominatedby = zeros(Int, numrow)
    frontier = []

    for i in 1:numrow
        for j in i+1:numrow
            row1 = df[i, :]
            row2 = df[j, :]
            row1count = 0
            row2count = 0
            for k in 1:numcol
                if row1[k] < row2[k]
                    row1count += 1
                elseif row1[k] > row2[k]
                    row2count += 1
                end
            end

            if row1count > 0 && row2count == 0
                dominatedby[j] += 1
                push!(dominates[i], j)
            elseif row2count > 0 && row1count == 0
                dominatedby[i] += 1
                push!(dominates[j], i)
            end
        end

        if dominatedby[i] == 0
            push!(frontier, i)
        end
    end

    rank = 1
    while (true)
        nextFrontier = []

        while !isempty(frontier)
            elem = popfirst!(frontier)
            endrank[elem] = rank
            for dominated in dominates[elem]
                dominatedby[dominated] -= 1
                if dominatedby[dominated] == 0
                    push!(nextFrontier, dominated)
                end
            end
        end

        if isempty(nextFrontier)
            break
        end

        frontier = nextFrontier
        rank += 1
    end
    
    DataFrame(RankColumn => endrank)
end

end