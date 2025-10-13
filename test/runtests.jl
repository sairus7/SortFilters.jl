using SortFilters
using Test
using Random
using SafeTestsets

# @test [] == detect_ambiguities(Base, Core)

@testset "SortFilters" begin

include("test_QuantileTracker.jl")
# @safetestset "movsort" begin include("test_movsort.jl") end
# @safetestset "movsort_stateful" begin include("test_movsort_stateful.jl") end

end
