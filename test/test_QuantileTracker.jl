x = rand(4)
xi = maximum(x)
qt = QuantileTracker(x, 1)
@test get(qt) === xi
add!(qt, 2.0)
@test get(qt) === 2.0

function test_consistency(::Type{T}; window_sizes, quantiles, n) where T
    inits = [rand(T, window_size) for window_size in window_sizes]
    # @show inits
    qts = [QuantileTracker(init, quantile) for quantile in quantiles, init in inits]
    msfs = [MovSortFilter(init) for init in inits]
    mismatches = 0
    for i in eachindex(window_sizes)
        for _ in 1:n
            value = rand(T)
            add!(msfs[i], value)
            for j in eachindex(quantiles)
                add!(qts[j, i], value)
                match = isapprox(get(msfs[i], quantiles[j]), get(qts[j, i]), rtol=eps(1.0))
                mismatches += !match
            end
        end
    end
    @test mismatches == 0
end

# TODO test NaN
test_consistency(
    Float64;
    window_sizes = vcat(1:6, rand(1:1000, 6)),
    quantiles = vcat(0, .5, 1, rand(6)),
    n = 500
)
