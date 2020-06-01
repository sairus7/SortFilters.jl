module SortFilters

export movsort!, movsort,
    MovSortFilter,
    reset!,
    add!, run, run!,
    get, get_max, get_min, get_median

##

const Evt{T} = NamedTuple{(:pos, :val),Tuple{Int,T}}
Evt(pos::Int, val::T) where T = (pos = pos, val = val)
tuplen(::NTuple{N, Any}) where {N} = N


# convert to rational type of p
movsort(in::AbstractVector{T}, window::Int, p::Union{Real, Tuple{Vararg{Real}}}) where T =
    movsort!(similar(in, typeof(p)), in, window, p)

# preserve type of input
movsort(in::AbstractVector{T}, window::Int, p::Integer) where T =
    movsort!(similar(in), in, window, p)

# preserve type of input
movsort(in::AbstractVector{T}, window::Int, p::Tuple{Vararg{Integer}}) where T =
    movsort!(similar(Array{NTuple{tuplen(p),T}}, axes(in)), in, window, p)

# sort vector in moving window
function movsort!(out, in::AbstractVector{T}, window::Int, p) where T
    len = length(p)
    buf = Vector{Evt{T}}(undef, window)
    fill!(buf, Evt(0, in[1]))

    for i in 1:length(in)
        x = in[i]

        # find insert position
        iold = 1
        inew = x > buf[1].val ? 1 : 0
        for i = 2 : window
            if buf[i].pos < buf[iold].pos
                iold = i
            end
            if x > buf[i].val
                inew = i
            end
        end

        # sort new element by memory shift and insert
        if inew > iold
            copyto!(buf, iold, buf, iold + 1, inew - iold)
        elseif iold > inew
            inew += 1
            copyto!(buf, inew + 1, buf, inew, iold - inew)
        end
        buf[inew] = Evt(i, x)

        out[i] = _quantile(buf, p)
    end
    out
end

mutable struct MovSortFilter{T}
    need_restart::Bool
    counter::Int
    buf::Vector{Evt{T}}
    MovSortFilter{T}(window::Int) where {T} =
        new{T}(true, 0, Vector{Evt{T}}(undef, window))
end

function reset!(state::MovSortFilter, window::Int)
    resize!(buf, window)
    reset!(state)
end

function reset!(state::MovSortFilter)
    state.need_restart = true
    state
end

function _restart!(state::MovSortFilter{T}, x0::T) where T
    fill!(state.buf, Evt(0, x0))
    state.counter = 1 # after first point
    state.need_restart = false
    state
end

# add one point
function add!(state::MovSortFilter{T}, x::T) where T
    if state.need_restart
        _restart!(state, x)
        return state
    end
    iold = 1
    buf = state.buf
    inew = x > buf[1].val ? 1 : 0
    for i = 2 : length(buf)
        if buf[i].pos < buf[iold].pos
            iold = i
        end
        if x > buf[i].val
            inew = i
        end
    end

    counter = state.counter
    # sort new element by memory shift and insert
    if inew > iold
        copyto!(buf, iold, buf, iold + 1, inew - iold)
    elseif iold > inew
        inew += 1
        copyto!(buf, inew + 1, buf, inew, iold - inew)
    end
    buf[inew] = Evt(counter, x)
    counter += 1
    state.counter = counter
    return state
end

# preserve original type
@inline get_max(state::MovSortFilter) = state.buf[end].val
@inline get_min(state::MovSortFilter) = state.buf[1].val

@inline get_median(state::MovSortFilter) = _median(state.buf)
@inline Base.get(state::MovSortFilter, p = 0.5) = _quantile(state.buf, p)
@inline function Base.get(state::MovSortFilter, p::Union{Real, Tuple{Vararg{Real}}})
    map(x->_quantile(state.buf, x), p)
end

# convert to rational type of p
function Base.run(state::MovSortFilter{T}, x::AbstractVector{T}, p::Union{Real, Tuple{Vararg{Real}}}) where T
    y = similar(x, typeof(p))
    run!(y, state, x, p)
end
# preserve type of input & simplified quantiles as nearest index
function Base.run(state::MovSortFilter{T}, x::AbstractVector{T}, p::Integer) where T
    y = similar(x)
    run!(y, state, x, p)
end
# preserve type of input & simplified quantiles as nearest index
function Base.run(state::MovSortFilter{T}, x::AbstractVector{T}, p::Tuple{Vararg{Integer}}) where T
    y = similar(Array{NTuple{tuplen(p),T}}, axes(x))
    run!(y, state, x, p)
end

function run!(
    y::AbstractVector, # lets not over-specialize ::AbstractVector{Union{T, Real, Tuple{Vararg{Real}}}}, # syntax highlight fails...
    state::MovSortFilter{T},
    x::AbstractVector{T},
    p # lets not over-specialize ::Union{Real, Tuple{Vararg{Real}}}
) where {T}
    length(x) == length(y) || throw(DimensionMismatch("arrays must have equal length"))
    for i in eachindex(x)
        y[i] = get(add!(state, x[i]), p)
    end
    y
end

## functions adopted from Statistics.jl, assumes `v` sorted

# converts type to rational
function _median(v::Vector{Evt{T}}) where T
    n = length(v)
    mid = div(n + 1, 2)
    if isodd(n)
        return v[mid].val |> Float64
    else
        return v[mid].val / 2 + v[mid+1].val / 2
    end
end

# converts type to rational
@inline function _quantile(v::Vector{Evt{T}}, p::Real) where T
    0 <= p <= 1 || throw(ArgumentError("input probability out of [0,1] range"))
    # require_one_based_indexing(v)

    lv = length(v)
    f0 = (lv - 1)*p # 0-based interpolated index
    t0 = trunc(f0)
    h  = f0 - t0
    i  = trunc(Int,t0) + 1

    a = v[i].val
    b = v[i + (h > 0)].val
    if isfinite(a) && isfinite(b)
        return a + h*(b-a)
    else
        return (1-h)*a + h*b
    end
end
# preserves original type when p is Int in 0..100 - just return nearest index
@inline function _quantile(v::Vector{Evt{T}}, p::Integer) where T
    0 <= p <= 100 || throw(ArgumentError("input probability out of [0,100] percent range"))
    # iq = (length(v) - 1) * p ÷ 100 + 1 # lower index
    iq = 1 + round(Int, (length(v) - 1) * p / 100, RoundNearestTiesUp) # nearest index
    return v[iq].val
end
@inline function _quantile(v::Vector{Evt{T}}, p) where T
    map(pi->_quantile(v, pi), p)
end

end # module
