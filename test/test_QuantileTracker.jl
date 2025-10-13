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
        last_nan = window_sizes[i] + 1
        for _ in 1:n
            value = rand() < .001 ? T(NaN) : rand(T)
            last_nan = isnan(value) ? last_nan + 1 : 0
            add!(msfs[i], value)
            for j in eachindex(quantiles)
                add!(qts[j, i], value)
                msf_val = get(msfs[i], quantiles[j])
                qt_val = get(qts[j, i])
                # MovSortFilter has unspecified behavior in the presence of NaN, so don't compare.
                # QuantileTracker sorts NaN at end line Base.sort.
                match = last_nan <= window_sizes[i] || isapprox(msf_val, qt_val, rtol=eps(1.0))
                mismatches += !match
                match || println("Mismatch at i=$i, j=$j: value=$value, quantile=$(quantiles[j]), window_size=$(window_sizes[i]), msf=$msf_val, qt=$qt_val")
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
