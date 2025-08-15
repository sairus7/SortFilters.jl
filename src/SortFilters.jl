module SortFilters

export movsort!, movsort,
    MovSortFilter, QuantileTracker,
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

function MovSortFilter(x::AbstractVector)
    msf = MovSortFilter{eltype(x)}(length(x))
    for value in x
        add!(msf, value)
    end
    msf
end

struct QuantileTracker{T, I}
    index::I # The index that this tracker is tracking, e.g. the 73rd element (discrete index, not a quantile)
    window_head::Base.RefValue{I} # Head of the circular queue `window`
    heaps::Memory{Tuple{T, I}} # Two heaps, a min heap below and a max heap above the quantile. Binary index trees.
    window::Memory{I} # Pointer to the current position in the heap of each value, stored in insertion order, circular fifo queue
    function QuantileTracker{T, I}(data::AbstractVector{T}; quantile=nothing, index=round(I, (length(data)-1)*quantile+1)) where {T, I}
        Base.require_one_based_indexing(data)
        heaps = Memory{Tuple{T, I}}(undef, 2length(data)+3) # Extra room to pad with typemax to avoid boundschecking

        checkbounds(data, index)

        # Copy data into the heap
        for i in eachindex(data)
            heaps[i+index+1] = (data[i], i) # TODO pack into a single UInt
        end

        # Heapify data
        sort!(view(heaps,index+2:index+length(data)+1), by=first)

        # Pad with typemax and typemin to avoid boundschecking (bubbling will never propagate through these values)
        for i in 1:index+1
            heaps[i] = (typemin(T), zero(I))
        end
        for i in index+length(data)+2:lastindex(heaps)
            heaps[i] = (typemax(T), zero(I))
        end

        # Populate window
        window = Memory{I}(undef, length(data))
        window_head = 1
        for i in index+2:index+length(data)+1
            window[heaps[i][2]] = i
        end

        new{T, I}(2index+1, Ref(window_head), heaps, window)
    end
end
QuantileTracker(data::AbstractVector{T}; kw...) where T = QuantileTracker{T, Int}(data; kw...)

# binary index tree indexing arithmetic
bit_parent(i) = i >> 1
bit_left_child(i) = i << 1
bit_right_child(i) = (i << 1) + one(i)

Base.get(qt::QuantileTracker) = qt.heaps[qt.index][1]

function _setindex!(qt, x, j)
    qt.heaps[j] = x
    qt.window[x[2]] = j
    nothing
end
# cmp(child, parent) is in order
function bubble_parent(get_parent, cmp, qt, value, j, jp, x) # <
    while true
        _setindex!(qt, x, j)
        j = jp
        j == qt.index && return j
        jp = get_parent(j)
        x = qt.heaps[jp]
        cmp(value, x[1]) && return j
    end
end
# cmp(child, parent) is in order
function bubble_child(get_left_child, cmp, qt, value, j)
    while true
        j_lc = get_left_child(j)
        x_l = qt.heaps[j_lc] # This being inbounds depends on the buffer space filled with typemax/typemin.
        j_rc = j_lc + one(j_lc)
        x_r = qt.heaps[j_rc]
        jc, x = cmp(x_l[1], x_r[1]) ? (j_rc, x_r) : (j_lc, x_l) # Only compare to the more parent like child
        cmp(x[1], value) && return j # the only end condition is the bubble reaching its appropraite location
        _setindex!(qt, x, j)
        j = jc
    end
end

function add!(qt::QuantileTracker{T}, value::T) where T
    # Legend:
    # i is qt.index
    # j is the index we are currently looking at
    # x is an element of interest in the heap (value, index)

    window_head = qt.window_head[]
    qt.window_head[] = mod(window_head + 1, eachindex(qt.window)) # Increment the head of the circular queue
    j = qt.window[window_head]
    # Replace oldest value in heaps with the new value, without relocating it
    # qt.heaps[j] = (value, j) (but don't actually do it, for performance)

    # Bubble the heap as needed
    i = qt.index
    if j > i # we're in the hi heap
        jp = bit_parent(j-i+1)+i-1
        x = qt.heaps[jp]
        if !(value >= x[1]) # Value is less than parent (out of order)
            j = bubble_parent(j -> bit_parent(j-i+1)+i-1, >=, qt, value, j, jp, x)
            if j == qt.index # we percolated all the way to the middle
                j = bubble_child(j -> i-bit_left_child(i-j+1), <=, qt, value, j)
            end
        else
            j = bubble_child(j -> bit_left_child(j-i+1)+i-1, >=, qt, value, j)
        end
    elseif j < i # we're in the lo heap
        jp = i-bit_parent(i-j+1)+1
        x = qt.heaps[jp]
        if !(value <= x[1]) # Value is less than parent (out of order)
            j = bubble_parent(j -> i-bit_parent(i-j+1)+1, <=, qt, value, j, jp, x)
            if j == qt.index # we percolated all the way to the middle
                j = bubble_child(j -> bit_left_child(j-i+1)+i-1, >=, qt, value, j)
            end
        else
            j = bubble_child(j -> i-bit_left_child(i-j+1), <=, qt, value, j)
        end
    else # we're at the middle
        j = bubble_child(j -> bit_left_child(j-i+1)+i-1, >=, qt, value, j) # try to bubble high
        if j == i # We didn't bubble high. Try to bubble low
            j = bubble_child(j -> i-bit_left_child(i-j+1), <=, qt, value, j)
        end
    end

    # Put the new value where the bubble ended up
    _setindex!(qt, (value, window_head), j)
end

end # module
