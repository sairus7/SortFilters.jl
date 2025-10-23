# to run the script, you should have these packages installed
using Plots
using StatsPlots
using BenchmarkTools
using FastRunningMedian
using Statistics
using SortFilters

y = randn(1000000)
windows = [5, 15, 51, 151, 515, 1515]
t_sf = zeros(size(windows))
t_qt = zeros(size(windows))
t_rm = zeros(size(windows))

function movquant(x::AbstractVector, window::Int, quantile::Real)
    qt = QuantileTracker(fill(x[1], window), quantile)
    mov_q = map(x) do xi
        add!(qt, Float64(xi))
        yi = get(qt)
    end
end

for i in eachindex(windows)
    w = windows[i]

    t = @benchmark SortFilters.movsort($y, $w, 0.5) samples=5
    t_sf[i] = Statistics.median(t).time / 10^6

    t = @benchmark movquant($y, $w, 0.5) samples=5
    t_qt[i] = Statistics.median(t).time / 10^6

    t = @benchmark running_median($y, $w) samples=5
    t_rm[i] = Statistics.median(t).time / 10^6
end

ctg = repeat([
    "1 - SortFilters . movsort",
    "2 - QuantileTracker",
    "3 - FastRunningMedian"], inner = length(windows))

# bars = log10.(hcat(t_sf, t_qt, t_rm)) # log scale
bars = hcat(t_sf, t_qt, t_rm)
xnames = repeat(1:length(windows), outer = size(bars,2))

groupedbar(xnames, bars, group = ctg, xlabel = "Window lengths", ylabel = "Time, ms",
        title = "Moving median performance, array length = 1000000.", bar_width = 0.7,
        lw = 0, framestyle = :box, legend=:topleft, ylims = (0, 500))

xticks!(1:length(windows), string.(windows))
# yticks!(collect(0:1:3), "10^{" .* string.(collect(0:1:3)) .* "}") # log scale

png("plot_qt.png")
