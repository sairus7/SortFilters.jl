using SortFilters
using Test
using Random
using Statistics


len = 100
Random.seed!(0)
T = Int
x = rand(T(-10):T(10), len)
window = 9

mx_ref = similar(x)
mn_ref = similar(x)
med_ref = similar(x, Float64)
p20_ref = similar(x, Float64)
p25_p75_ref = similar(x, NTuple{2,Float64})

for i = 1:len
    chunk = x[max(1, i-window+1) : i]
    mx_ref[i] = maximum(chunk)
    mn_ref[i] = minimum(chunk)
    med_ref[i] = median(chunk)
    p20_ref[i] = quantile(chunk, 0.20)
    p25_p75_ref[i] = quantile(chunk, (0.25, 0.75))
end


@testset "movsort" begin

mx = movsort(x, window, 1.0)
mx_ = movsort(x, window, 100)

mn = movsort(x, window, 0.0)
mn_ = movsort(x, window, 0)

med = movsort(x, window, 0.5)
med_ = movsort(x, window, 50)

p25_p75 = movsort(x, window, (0.25, 0.75))
p25_p75_ = movsort(x, window, (25, 75))

@test mx == mx_ == mx_ref
@test mn == mn_ == mn_ref
@test med[window:end] == med_[window:end] == med_ref[window:end]
@test p25_p75[window:end] == p25_p75_[window:end] == p25_p75_ref[window:end]

p20 = movsort(x, window, 0.20)
@test p20[window:end] == p20_ref[window:end] # should match nearest index quantiles

end


@testset "movsort!" begin

mx 		 = similar(mx_ref, Float64)
mx_ 	 = similar(mx_ref)
mn 		 = similar(mn_ref, Float64)
mn_ 	 = similar(mn_ref)
med 	 = similar(med_ref, Float64)
med_ 	 = similar(med_ref)
p25_p75  = similar(p25_p75_ref, NTuple{2,Float64})
p25_p75_ = similar(p25_p75_ref)
p20      = similar(p20_ref, Float64)

movsort!(mx, x, window, 1.0)
movsort!(mx_, x, window, 100)

movsort!(mn, x, window, 0.0)
movsort!(mn_, x, window, 0)

movsort!(med, x, window, 0.5)
movsort!(med_, x, window, 50)

movsort!(p25_p75, x, window, (0.25, 0.75))
movsort!(p25_p75_, x, window, (25, 75))

@test mx == mx_ == mx_ref
@test mn == mn_ == mn_ref
@test med[window:end] == med_[window:end] == med_ref[window:end]
@test p25_p75[window:end] == p25_p75_[window:end] == p25_p75_ref[window:end]

movsort!(p20, x, window, 0.20)
@test p20[window:end] == p20_ref[window:end] # should not match nearest index quantiles

end
