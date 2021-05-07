# SortFilters

Moving quantiles implemented as fast moving window sort algorithm.
Implemented both as functions over a moving window, and stateful filter objects.

Setting an appropriate probability level, you can get: moving median, minimum, maximum, quartiles and so on.

## Installation
```julia
pkg> add SortFilters
```

## Comparison with other packages
There are two other Julia packages with overlapping functionality for moving quantiles functions:
- [RollingFunctions.jl](https://github.com/JeffreySarnoff/RollingFunctions.jl)
- [ImageFiltering.jl](https://github.com/JuliaImages/ImageFiltering.jl)

Compared to these packages, [SortFilters.jl](https://github.com/sairus7/SortFilters.jl) provides significant speed-up (3x-20x faster), depending on window size (benchmark available at [examples/benchmark.jl](https://github.com/sairus7/SortFilters.jl/blob/master/examples/benchmark.jl)):

![plot](https://user-images.githubusercontent.com/20798349/83455606-e42d1280-a466-11ea-8b48-0e375a7dcaa0.png)

Also [SortFilters.jl](https://github.com/sairus7/SortFilters.jl) provides stateful filter objects, allowing you to process a signal of indefinite length in RAM-friendly chunks, similar to [DSP.jl](https://juliadsp.github.io/DSP.jl/stable/filters/#stateful-filter-objects-1).

However, if all you need is moving maximum, minimum, or range, then you should use another package, [MaxMinFilters.jl](https://github.com/sairus7/MaxMinFilters.jl), which implement much faster algorithms than moving sort.

## Examples
[examples/example.jl](https://github.com/sairus7/SortFilters.jl/blob/master/examples/example.jl):
```julia
using Plots
using SortFilters
using Random

Random.seed!(5)
len = 200
x = trunc.(Int, 100*randn(len))
x = cumsum(x)
w = 20

# when p represents probability rational number in the range of [0..1]
probability = 0.5
# output is converted to rational type of p
moving_median = movsort(x, w, probability)

# when p represents integer percent number in the range of [0..100]
percent = 50
# output preserves input type and is calculated as simplified quantiles
# from input elements at nearest index
moving_median_int = movsort(x, w, percent)

# also you can get several quantiles at once:
p = (0.25, 0.75)
moving_quartiles = movsort(x, w, p)
# and then calculate interquantile range
moving_iqr = map( x->x[2]-x[1], moving_quartiles)
moving_q25 = map( x->x[1], moving_quartiles)
moving_q75 = map( x->x[2], moving_quartiles)

# notice that quantiles are calculated with a delay of window/2
# start boundary condition is set as a repeated first point
plot(x, label = "x", legend = :topleft)
plot!(moving_median, label = "exact moving median")
plot!(moving_median_int, label = "simplified moving median")
plot!(moving_q25, label = "25th quartile")
plot!(moving_q75, label = "75th quartile")

#png("plot1.png")
```
![plot1](https://user-images.githubusercontent.com/20798349/83455610-e68f6c80-a466-11ea-8fb3-7b0c18aa92d0.png)
