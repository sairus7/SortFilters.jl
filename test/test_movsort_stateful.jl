using SortFilters
using Test
using Random
using Statistics

#@testset "Stateful sort filters" begin

len = 100
Random.seed!(1)
T = Int
v = rand(0:100, len)
window = 5

mx_ref = similar(v)
mn_ref = similar(v)
med_ref = similar(v, Float64)
p03_ref = similar(v, Float64)
p02_p07_ref = similar(v, NTuple{2,Float64})

for i = 1:len
    chunk = v[max(1, i-window+1) : i]
    mx_ref[i] = maximum(chunk)
    mn_ref[i] = minimum(chunk)
    med_ref[i] = median(chunk)
    p03_ref[i] = quantile(chunk, 0.3)
    p02_p07_ref[i] = quantile(chunk, (0.2, 0.7))
end

#Juno.@enter median(v[1:2])
#Juno.@enter quantile(v[1:2], 0.5)
#@testset "MovSortFilter" begin
f = MovSortFilter{T}(window)
mn = run(f, v, 0)
reset!(f)
mx = run(f, v, 100)
reset!(f)
med = run(f, v, 0.5)
@test mn == mn_ref
@test mx == mx_ref
@test med[window:end] == med_ref[window:end] # it has differnt behaviour at start bound

# reset state
reset!(f)
med = run(f, v[1:len÷2], 0.5)
@test med[window:end] == med_ref[window:len÷2] # it has differnt behaviour at start bound
# maintain state
med = run(f, v[len÷2+1:end], 0.5)
@test med == med_ref[len÷2+1:end]


#end

#end



using SystemBenchmark

using JuliaDB

data_table = loadtable("C:\\Users\\gvg\\Downloads\\data.csv", delim = ',', header_exists = true)

data_table

using Plots

select(data_table, )
t = filter(row -> row.lang == "julia", data_table)
plot(t[Symbol(size(B)])

t = filter(row -> row.lang == "julia", data_table)

t[1]

t = filter(row -> getfield(t, Symbol("size(B)")) == "julia", data_table)

getfield


data = rand(1:100, 100)

using VegaLite

@vlplot(
[

]
    :line,
    x=:value,
    y=:value,
    row=:coordinate,
    column=:coordinate,
)

x = 1:100

y = rand(100)
[@vlplot(:line, x, y,); @vlplot(:point, x, y,)] |> display
