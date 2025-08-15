x = rand(4)
i = 4
xi = maximum(x)
qt = QuantileTracker(x; index=i)
@test get(qt) === xi
add!(qt, 2.0)
@test get(qt) === 2.0

function test_consistency(::Type{T}; window_sizes, quantiles, n) where T
    inits = [rand(T, window_size) for window_size in window_sizes]
    # @show inits
    qts = [QuantileTracker(init; quantile) for quantile in quantiles, init in inits]
    msfs = [MovSortFilter(init) for init in inits]
    mismatches = 0
    for i in eachindex(window_sizes)
        for _ in 1:n
            value = rand(T)
            add!(msfs[i], value)
            for j in eachindex(quantiles)
                add!(qts[j, i], value)
                q = window_sizes[i] == 1 ? quantiles[j] : (round(Int, quantiles[j] * (window_sizes[i]-1) + 1) - 1) / (window_sizes[i]-1)
                match = get(msfs[i], q) ≈ get(qts[j, i])
                mismatches += !match
                # @show qts[j, i] get(qts[j, i]) get(msfs[i], q)
                # get(msfs[i], q) !== get(qts[j, i]) && println("ERROR: ", (window_sizes[i],quantiles[j],q))
                # @test get(msfs[i], quantiles[j]) === get(qts[j, i]) Uncomment to fail eagerly
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
