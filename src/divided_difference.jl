###############
# API methods #
###############

"""
    div_diff(f, x::Vector; ill_test=true, kwargs...)
    div_diff(f, x::Tuple; kwargs...)
    div_diff(f, x...; kwargs...)

Return the divided difference `f[x_0,x_1,...,x_n]`, assuming `f` is called as `f(x)`. 
A ill-conditioned test for the matrix function computation is set by default. 
You can set `ill_test=false` to obtain a quik result. 
"""

@inline function div_diff(f::F, x::Vector{R}; 
                          ill_test=true, 
                          kwargs...) where {F,R<:Real}
    sort!(x)
    fd = mat_fun(f, FiniteDual(x); ill_test, kwargs...)
    return extract_DD(fd)
end
@inline div_diff(f::F, x::Tuple; kwargs...) where {F} = div_diff(f, collect(x); kwargs...)
@inline div_diff(f::F, x...; kwargs...) where {F} = div_diff(f, tuple(x...); kwargs...)

"""
    div_diff(f!, y::AbstractArray, x::Vector{Real}; kwargs...)

Return the divided difference `f![x_0,x_1,...,x_n]`, assuming `f!` is called as `f!(y, x)`, 
where `y` is an array and f!(x_0) is stored in `y`.
"""
@inline function div_diff(f!::F, y::AbstractArray{Y}, 
                          x::Vector{R}; kwargs...) where {F,Y,R<:Real}
    sort!(x)
    require_one_based_indexing(y)
    yfd = similar(y, FiniteDual{Y,length(x)})
    mat_fun!(f!, yfd, FiniteDual(x); kwargs...)
    v0(x) = value(x, 1)
    map!(v0, y, yfd)
    return extract_DD(yfd)
end

"""
    div_diff!(result::AbstractArray, f, x::Vector{Real})

Compute the divided difference `f[x_0,x_1,...,x_n]` and store the results in `result`, 
assuming `f` is called as `f(x)` and returns array.
"""
@inline function div_diff!(result::AbstractArray, f::F,
                           x::Vector{R}; kwargs...) where {F,R<:Real}
    sort!(x)
    require_one_based_indexing(result)
    yfd = mat_fun(f, FiniteDual(x); kwargs...)
    result = extract_DD!(result, yfd)
    return result
end

"""
    div_diff!(result::AbstractArray, f!, y::AbstractArray, x::Vector{Real})

Compute the divided difference `f![x_0,x_1,...,x_n]` and store the results in `result`, 
assuming `f!` is called as `f!(y, x)`, where `y` is an array and f!(x_0) is stored in `y`.
"""
@inline function div_diff!(result::AbstractArray, f!::F,
                           y::AbstractArray{Y}, 
                           x::Vector{R}; kwargs...) where {F,Y,R<:Real}
    sort!(x)
    require_one_based_indexing(result, y)
    yfd = similar(y, FiniteDual{Y,length(x)})
    mat_fun!(f!, yfd, FiniteDual(x); kwargs...)
    v0(x) = value(x, 1)
    map!(v0, y, yfd)
    result = extract_DD!(result, yfd)
    return result
end

#####################
# result extraction #
#####################

# non-mutating #
#--------------#

@inline extract_DD(y::Real)          = zero(y)
@inline extract_DD(y::Complex)       = zero(y)
@inline extract_DD(y::FiniteDual)    = y.DD_table[1, length(y)]
@inline extract_DD(y::AbstractArray) = map(dd -> extract_DD(dd), y)

# mutating #
#----------#

extract_DD!(result::AbstractArray, y::AbstractArray) = map!(dd -> extract_DD(dd), result, y)
