module testpareto

include("../src/pareto.jl")

using Test
using .Pareto
using DataFrames

names = ["a"; "b"]
df1 = DataFrame([0 1; 1 2], names)
df2 = DataFrame([1 1; 1 1], names)
df3 = DataFrame([1 1; 1 2], names)
df4 = DataFrame([2 2; 1 1], names)
df5 = DataFrame([], [])

@test calcranks(df1) == DataFrame(Pareto.RankColumn => [1,2])
@test calcranks(df2) == DataFrame(Pareto.RankColumn => [1,1])
@test calcranks(df3) == DataFrame(Pareto.RankColumn => [1,2])
@test calcranks(df4) == DataFrame(Pareto.RankColumn => [2,1])
@test calcranks(df5) == DataFrame(Pareto.RankColumn => [])

end