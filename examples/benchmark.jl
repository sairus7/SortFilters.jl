# to run the script, you should have these packages installed
using Plots
using StatsPlots
using BenchmarkTools
using Statistics
using Indicators
using RollingFunctions
using ImageFiltering
using SortFilters
using MaxMinFilters

y = randn(1000000)
windows = [5, 15, 51, 151]
t_mf1 = zeros(size(windows))
t_mf2 = zeros(size(windows))
t_img = zeros(size(windows))
t_rol = zeros(size(windows))
t_ind = zeros(size(windows))

w = windows[4]

@time out_ = MaxMinFilters.movmax(y, w)
@time out1 = SortFilters.movsort(y, w, 1)

y = randn(200)
@time out1 = SortFilters.movsort(y, w, (100, 50))

@time out2 = RollingFunctions.rollmedian(y, w)

@time out3 = ImageFiltering.mapwindow(median!, y, w)
n = 200
plot(out1[1:n])
plot!(out2[1:n])
plot!(out3[1:n])

# test for moving maximum

for i in eachindex(windows)
    w = windows[i]

    t = @benchmark MaxMinFilters.movmax($y, $w) samples=3
    t_mf1[i] = median(t).time / 10^6

    t = @benchmark MaxMinFilters.movmaxmin($y, $w) samples=3
    t_mf2[i] = median(t).time / 10^6

    t = @benchmark ImageFiltering.mapwindow($extrema, $y, $w) samples=3
    t_img[i] = median(t).time / 10^6

    t = @benchmark RollingFunctions.rollmax($y, $w) samples=3
    t_rol[i] = median(t).time / 10^6

    t = @benchmark Indicators.runmax($y, n = $w, cumulative = $false) samples=3
    t_ind[i] = median(t).time / 10^6
end

ctg = repeat(["1 - MaxMinFilters . movmax",
              "2 - MaxMinFilters . movmaxmin",
              "3 - ImageFiltering . mapwindow extrema",
              "4 - RollingFunctions . rollmax",
              "5 - Indicators . runmax",
              ], inner = length(windows))

#bars = log10.(hcat(t_mf1, t_mf2, t_img, t_rol, t_ind)) # log scale
bars = hcat(t_mf1, t_mf2, t_img, t_rol, t_ind)
xnames = repeat(collect(1:length(windows)), outer = size(bars,2))

groupedbar(xnames, bars, group = ctg, xlabel = "Window lengths", ylabel = "Time, ms",
        title = "Moving maximum performance, array length = 1000000.", bar_width = 0.7,
        lw = 0, framestyle = :box, legend=:topleft, ylims = (0, 530))

xticks!(collect(1:4), string.(windows))
#yticks!(collect(0:1:3), "10^{" .* string.(collect(0:1:3)) .* "}") # log scale

#png("plot.png")
