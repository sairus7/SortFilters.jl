# to run the script, you should have these packages installed
using Plots
using StatsPlots
using BenchmarkTools
using Statistics
using RollingFunctions
using ImageFiltering
using SortFilters

y = randn(1000000)
windows = [5, 15, 51, 151]
t_srt = zeros(size(windows))
t_img = zeros(size(windows))
t_rol = zeros(size(windows))

for i in eachindex(windows)
    w = windows[i]

    t = @benchmark SortFilters.movsort($y, $w, 0.5) samples=3
    t_srt[i] = median(t).time / 10^6

    t = @benchmark ImageFiltering.mapwindow($median, $y, $w) samples=3
    t_img[i] = median(t).time / 10^6

    t = @benchmark RollingFunctions.rollmedian($y, $w) samples=3
    t_rol[i] = median(t).time / 10^6
end

ctg = repeat(["1 - SortFilters . movsort",
              "2 - ImageFiltering . mapwindow median",
              "3 - RollingFunctions . rollmedian",
              ], inner = length(windows))

#bars = log10.(hcat(t_mfl, t_srt, t_img, t_rol, t_ind)) # log scale
bars = hcat(t_srt, t_img, t_rol)
xnames = repeat(1:length(windows), outer = size(bars,2))

groupedbar(xnames, bars, group = ctg, xlabel = "Window lengths", ylabel = "Time, ms",
        title = "Moving median performance, array length = 1000000.", bar_width = 0.7,
        lw = 0, framestyle = :box, legend=:topleft, ylims = (0, 5200))

xticks!(1:4, string.(windows))
#yticks!(collect(0:1:3), "10^{" .* string.(collect(0:1:3)) .* "}") # log scale

png("plot.png")
