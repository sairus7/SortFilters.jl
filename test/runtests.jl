using SortFilters
using Test
using Random

@test [] == detect_ambiguities(Base, Core)

tests = [
         "movsort_stateful",
        ]

if length(ARGS) > 0
    tests = ARGS
end

@testset "SortFilters" begin

for t in tests
    fp = joinpath(dirname(@__FILE__), "test_$t.jl")
    println("$fp ...")
    include(fp)
end

end

#=
using Statistics
v = [1,2,3,4,5,6,7,8,9,10,11]
Juno.@enter quantile!(v, 0.5; sorted=true)
Juno.@enter quantile!(v, (0.5, 0.7); sorted=true)
=#
