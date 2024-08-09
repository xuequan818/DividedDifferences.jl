###############
# API methods #
###############

"""
    divided_difference(f, x::Vector)
    divided_difference(f, x::Tuple)
    divided_difference(f, x...)

Return the divided difference `f[x_0,x_1,...,x_n]`, assuming `f` is called as `f(x)`.
"""
@inline function divided_difference(f::F, x::Vector{R}) where {F,R<:Real}
    FD = FiniteDual(x)
    return extract_DD(f(FD))
end
@inline divided_difference(f::F, x::Tuple) where {F} = divided_difference(f, collect(x))
@inline divided_difference(f::F, x...) where {F} = divided_difference(f, tuple(x...))

"""
    divided_difference(f!, y::AbstractArray, x::Vector{Real})

Return the divided difference `f![x_0,x_1,...,x_n]`, assuming `f!` is called as `f!(y, x)`, 
where `y` is an array and f!(x_0) is stored in `y`.
"""
@inline function divided_difference(f!::F, y::AbstractArray{Y}, 
                                    x::Vector{R}) where {F,Y,R<:Real}
    require_one_based_indexing(y)
    yfd = similar(y, FiniteDual{Y,length(x)})
    f!(yfd, FiniteDual(x))
    v0(x) = value(x, 1)
    map!(v0, y, yfd)
    return extract_DD(yfd)
end

"""
    divided_difference!(result::AbstractArray, f, x::Vector{Real})

Compute the divided difference `f[x_0,x_1,...,x_n]` and store the results in `result`, 
assuming `f` is called as `f(x)` and returns array.
"""
@inline function divided_difference!(result::AbstractArray, f::F,
                                     x::Vector{R}) where {F,R<:Real}
    require_one_based_indexing(result)
    yfd = f(FiniteDual(x))
    result = extract_DD!(result, yfd)
    return result
end

"""
    divided_difference!(result::AbstractArray, f!, y::AbstractArray, x::Vector{Real})

Compute the divided difference `f![x_0,x_1,...,x_n]` and store the results in `result`, 
assuming `f!` is called as `f!(y, x)`, where `y` is an array and f!(x_0) is stored in `y`.
"""
@inline function divided_difference!(result::AbstractArray, f!::F,
                                     y::AbstractArray{Y}, 
                                     x::Vector{R}) where {F,Y,R<:Real}
    require_one_based_indexing(result, y)
    yfd = similar(y, FiniteDual{Y,length(x)})
    f!(yfd, FiniteDual(x))
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
