using SortFilters
using Test
using Random
using Statistics


len = 100
Random.seed!(1)
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


@testset "reset filter state and run" begin

f = MovSortFilter{T}(window)
mx = run(f, x, 1.0)
reset!(f)
mx_ = run(f, x, 100)

reset!(f)
mn = run(f, x, 0.0)
reset!(f)
mn_ = run(f, x, 0)

reset!(f)
med = run(f, x, 0.5)
reset!(f)
med_ = run(f, x, 50)

reset!(f)
p25_p75 = run(f, x, (0.25, 0.75))
reset!(f)
p25_p75_ = run(f, x, (25, 75))

@test mx == mx_ == mx_ref
@test mn == mn_ == mn_ref
@test med[window:end] == med_[window:end] == med_ref[window:end]
@test p25_p75[window:end] == p25_p75_[window:end] == p25_p75_ref[window:end]

reset!(f)
p20 = run(f, x, 0.20)
@test p20[window:end] ≈ p20_ref[window:end] # should match nearest index quantiles

end


@testset "maintain filter state and run!" begin

mx 		 = similar(mx_ref, Float64)
mx_ 	 = similar(mx_ref)
mn 		 = similar(mn_ref, Float64)
mn_ 	 = similar(mn_ref)
med 	 = similar(med_ref, Float64)
med_ 	 = similar(med_ref)
p25_p75  = similar(p25_p75_ref, NTuple{2,Float64})
p25_p75_ = similar(p25_p75_ref)
p20      = similar(p20_ref, Float64)

part1 = x->view(x, 1 : len÷2)
part2 = x->view(x, len÷2+1 : len)
part1(x)

f = MovSortFilter{T}(window)
run!(part1(mx), f, part1(x), 1.0)
run!(part2(mx), f, part2(x), 1.0)
reset!(f)
run!(part1(mx_), f, part1(x), 100)
run!(part2(mx_), f, part2(x), 100)

reset!(f)
run!(part1(mn), f, part1(x), 0.0)
run!(part2(mn), f, part2(x), 0.0)
reset!(f)
run!(part1(mn_), f, part1(x), 0)
run!(part2(mn_), f, part2(x), 0)

reset!(f)
run!(part1(med), f, part1(x), 0.5)
run!(part2(med), f, part2(x), 0.5)
reset!(f)
run!(part1(med_), f, part1(x), 50)
run!(part2(med_), f, part2(x), 50)

reset!(f)
run!(part1(p25_p75), f, part1(x), (0.25, 0.75))
run!(part2(p25_p75), f, part2(x), (0.25, 0.75))
reset!(f)
run!(part1(p25_p75_), f, part1(x), (25, 75))
run!(part2(p25_p75_), f, part2(x), (25, 75))

@test mx == mx_ == mx_ref
@test mn == mn_ == mn_ref
@test med[window:end] == med_[window:end] == med_ref[window:end]
@test p25_p75[window:end] == p25_p75_[window:end] == p25_p75_ref[window:end]

reset!(f)
run!(part1(p20), f, part1(x), 0.20)
run!(part2(p20), f, part2(x), 0.20)
@test p20[window:end] ≈ p20_ref[window:end] # should match nearest index quantiles

end
