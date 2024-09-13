import Base: +, -, *, /, ^, >, <, ==, >=, <=

##############
# FiniteDual #
##############

"""
    DividedDifferences.can_finitedual(V::Type)

Determines whether the type V is allowed as the scalar type in a
FiniteDual. By default, only `<:Real` types are allowed.
"""

can_finitedual(::Type{<:Number}) = true
can_finitedual(::Type) = false

"""
    DividedDifferences.can_dd_table(DD_table::AbstractMatrix)

Test whether the DD_table is an upper triangular matrix.
"""

can_dd_table(A::AbstractMatrix) = istriu(A) || norm(A-UpperTriangular(A)) < sqrt(eps(eltype(A)))

struct FiniteDual{T,N} <: Number
    DD_table::SMatrix{N, N, T}
    function FiniteDual(DD_table::SMatrix{N,N,T}) where {T,N}
        can_finitedual(T) || throw_cannot_finitedual(T) 
        can_dd_table(DD_table) || throw_cannot_dd_table(DD_table) 
        new{T,N}(DD_table)
    end
end

##############
# Exceptions #
##############

@noinline function throw_cannot_finitedual(T::Type)
    throw(ArgumentError("Cannot create a finite dual over scalar type $T." *
                        " If the type behaves as a scalar, define DividedDifferences.can_finitedual(::Type{$T}) = true."))
end

@noinline function throw_cannot_dd_table(A::AbstractMatrix)
    error("Divided difference table should be upper triangular.")
end

################
# Constructors #
################

@inline function FiniteDual(x::AbstractMatrix{T}, ::Val{N}) where {T,N}
    N == 1 ? (can_finitedual(T) ? x : throw_cannot_finitedual(T)) : FiniteDual(SMatrix{N,N,T}(x))
end
@inline FiniteDual(x::AbstractMatrix) = FiniteDual(x, Val(size(x,1)))

# construct FiniteDual for the identity function f(x) = x
@inline FiniteDual(x::AbstractVector{T}) where {T} = FiniteDual(Bidiagonal{T}(x, ones(T,(length(x) - 1)), :U))
@inline function FiniteDual(x::Tuple)
    T = promote_type(typeof.(x)...)
    FiniteDual(collect(T,x))
end
@inline FiniteDual(x...) = FiniteDual(tuple(x...))

# construct FiniteDual for the constant function f(x) = c
@inline FiniteDual{N}(c::Number) where {N} = FiniteDual{N}(c*I)
@inline FiniteDual{N}(J::UniformScaling) where {N} = FiniteDual(SDiagonal{N}(J))

##############################
# Utility/Accessor Functions #
##############################

@inline Base.eltype(::FiniteDual{T,N}) where {T,N} = T
@inline Base.eltype(::Type{FiniteDual{T,N}}) where {T,N} = T

@inline Base.length(::FiniteDual{T,N}) where {T,N} = N
@inline Base.length(::Type{FiniteDual{T,N}}) where {T,N} = N

@inline table(x::FiniteDual) = x.DD_table

@inline Base.values(x::FiniteDual) = diag(table(x))
@inline function value(x::FiniteDual{T,N}, i::Integer) where {T,N}
    @assert 1 <= i <= N
    return values(x)[i]
end

########################
# Promotion/Conversion #
########################

function Base.promote_rule(::Type{FiniteDual{T,N}},
                           ::Type{FiniteDual{V,N}}) where {T,V,N}
    return FiniteDual{promote_type(T, V),N}
end

for R in (AbstractIrrational, Real, BigFloat, Bool)
    if isconcretetype(R)
        @eval begin
            Base.promote_rule(::Type{$R}, ::Type{FiniteDual{T,N}}) where {T,N} = FiniteDual{promote_type($R, T),N}
            Base.promote_rule(::Type{FiniteDual{T,N}}, ::Type{$R}) where {T,N} = FiniteDual{promote_type($R, T),N}
        end
    else
        @eval begin
            Base.promote_rule(::Type{R}, ::Type{FiniteDual{T,N}}) where {R<:$R,T,N} = FiniteDual{promote_type(R, T),N}
            Base.promote_rule(::Type{FiniteDual{T,N}}, ::Type{R}) where {T,N,R<:$R} = FiniteDual{promote_type(R, T),N}
        end
    end
end

@inline Base.convert(::Type{FiniteDual{T,N}}, x::FiniteDual{V,N}) where {T,V,N} = FiniteDual(T.(table(x)))
@inline Base.convert(::Type{FiniteDual{T,N}}, x::Number) where {T,N} = FiniteDual{N}(T(x))
@inline Base.convert(::Type{FD}, x::FD) where {FD<:FiniteDual} = x

#####################
# Generic Functions #
#####################

@inline Base.zero(x::FiniteDual) = zero(typeof(x))
@inline Base.zero(::Type{FiniteDual{T,N}}) where {T,N} = FiniteDual(zeros(SMatrix{N,N,T}))

@inline Base.one(x::FiniteDual) = one(typeof(x))
@inline Base.one(::Type{FiniteDual{T,N}}) where {T,N} = FiniteDual{N}(T(1))

# Predicates #
#------------#

@inline isconstant(x::FiniteDual) = allequal(values(x))

@inline Base.iszero(x::FiniteDual) = iszero(table(x))

for pred in [:>, :<, :(==), :>=, :<=]
    @eval @inline $pred(x::FiniteDual, y::Number) = 
                   build_pred_table(x, $pred.(values(x), y))
    @eval @inline $pred(x::Number, y::FiniteDual) = 
                   build_pred_table(y, $pred.(x, values(y)))
end

# The prediction function `f` acts on a matrix `A` could be computed by 
# solving the Sylvester equation for 2×2 block matrices. 
# Compute `f(A)` by the Schur decomposition `f(A) = Z*f(T)*Z'`.
# First, reorder and divide `T` into a 2×2 matrix of upper triangular blocks, 
# where one of the diagonal blocks has all eigenvalues meeting the prediction
# and the other does not. The corresponding diagonal parts of `F := f(T)` are 
# the unit and zero matrices. By the commutativity relation FT = TF, 
# there is the Sylvester equation: `T11*F12 - F12*T22 + (T12*F22-F11*T12) = 0`.
@inline function build_pred_table(x::FiniteDual{Tx,N}, pred_val) where {Tx,N}
    S = schur(table(x))
    block_size = pred_val[1] ? [sum(pred_val), N] : [N - sum(pred_val), N]
    fx = Array(eltype(S).(pred_val))
    if !(issorted(pred_val;rev=true) || issorted(pred_val;rev=false))
        S, block_size = reorder_schur(schur(table(x)), Int.(pred_val) .+ 1)
        sort!(fx; rev=true)
    end
    T, Z, Λ = S

    F = diagm(fx)
    i = 1:block_size[1]
    j = block_size[1]+1:block_size[2]
    Y = T[i, j] * F[j, j] - F[i, i] * T[i, j]

    # solve Tii*Fij - Fij*Tjj + Y = 0
    if length(i) > 1 || length(j) > 1
        F[i, j] = sylvester(T[i, i], -T[j, j], Y)
    else
        F[i, j] = Y ./ (T[i, i] - T[j, j])
    end

    return FiniteDual(Z * F * Z')
end

####################
# Matrix Functions #
####################

"""
    DividedDifferences.mat_fun(f, x::FiniteDual; ill_test::Bool=true, kwargs...)

Compute the matrix function `f(x.DD_table)`. By default, 
the calculation is first performed by Julia matrix operations 
and then tested whether the result is ill-conditioned, 
and if it is ill-conditioned, the calculation is performed by the schur parlett algorithm. 
You can set `ill_test=false` to skip the testing and recalculation process.
"""

@inline function mat_fun(f::Function, x::FiniteDual{T,N};
                         ill_test::Bool=true, 
                         kwargs...) where {T,N}
    F = f(x)
    issign(F) && return extract_sign(F)
    N == 2 || isconstant(x) && return F
    
    ill_test ? (
        try
            _mat_fun(f, F, x, Val(ill_test))
        catch e
            if isa(e, TypeError)
                f!(y, x) = copy!(y, f(x))
                _mat_fun!(f!, F, x, Val(ill_test))
            else
                throw(e)
            end
        end
    ) : return F
end

@inline function _mat_fun(f::Function, F, x::FiniteDual,
    ::Val{true}; kwargs...)
    if isill_mat_fun(F)
        F = FiniteDual(mat_fun_schur_parlett(f, table(x); kwargs...))
    end
    return F
end

"""
    DividedDifferences.mat_fun!(f!, F::AbstractArray{T}, x::FiniteDual; ill_test::Bool=true, kwargs...) where {T<:FiniteDual}

Compute the matrix function `f(x.DD_table)`, and store the results in `F`. 
See `DividedDifferences.mat_fun` for a description of possible keyword arguments.
"""

@inline function mat_fun!(f!::Function, F::AbstractArray,
                          x::FiniteDual{T,N}; ill_test::Bool=true,
                          kwargs...) where {T,N}
    f!(F, x)
    N == 2 || isconstant(x) && return F

    ill_test ? _mat_fun!(f!, F, x, Val(ill_test); kwargs...) : F
end

@inline function _mat_fun!(f!::Function, F::AbstractArray{T},
                           x::FiniteDual, ::Val{true};
                           kwargs...) where {T<:FiniteDual}
    select = isill_mat_fun(F)
    if any(isone, select)
        f!_select(y, x) = begin
            f!(y, x)
            y[select]
        end
        F_select = Array.(table.(F[select]))
        F[select] = FiniteDual.(mat_fun_schur_parlett!(f!_select, F_select,
            table(x); f!_T=eltype(T), f!_S=size(F), kwargs...))
    end
    return F
end

@inline function isill_mat_fun(F::FiniteDual{V}) where {V}
    ill_cond = sqrt(eps(real(float(oneunit(V)))))
    Tf = table(F)
    (any(isnan, Tf) || any(isinf, Tf) || cond(Tf) > ill_cond) ? true : false
end
@inline isill_mat_fun(F::Number) = false
@inline isill_mat_fun(F::AbstractArray) = map(x -> isill_mat_fun(x), F)

###################################
# General Mathematical Operations #
###################################

#################
# Special Cases #
#################

# (+)/(-) #
#---------#

for op in (:+, :-)
    @eval @inline function $op(x::FiniteDual, y::FiniteDual)
        @assert length(x) == length(y)
		FiniteDual($op(table(x),table(y)))
	end
    @eval @inline $op(x::FiniteDual{Tx,Nx}, y::Number) where {Tx,Nx} = FiniteDual($op(table(x), SDiagonal{Nx}(y * I)))
    @eval @inline $op(x::Number, y::FiniteDual{Ty,Ny}) where {Ty,Ny} = FiniteDual($op(SDiagonal{Ny}(x * I), table(y)))
end
@inline -(x::FiniteDual) = FiniteDual(-table(x))

# (*)/(/) #
#---------#

for op in (:*, :/)
    @eval @inline function $op(x::FiniteDual, y::FiniteDual)
        @assert length(x) == length(y)
        FiniteDual($op(table(x), table(y)))
    end
    @eval @inline $op(x::FiniteDual, y::Number) = FiniteDual($op(table(x), y))
end
@inline *(x::Number, y::FiniteDual) = FiniteDual(x * table(y))
@inline /(x::Number, y::FiniteDual) = FiniteDual(x * inv(table(y)))

# matrix functions #
#------------------#

for f in (:exp, :cis, :sqrt,
          :asin, :acos, :atan,
          :acsc, :asec, :acot,
          :asinh, :acosh, :atanh,
          :acsch, :asech, :acoth)
    @eval @inline Base.$f(x::FiniteDual) = FiniteDual($f(table(x)))
end

for f in (:log, :sin, :cos, :tan,
          :csc, :sec, :cot,
          :sinh, :cosh, :tanh,
          :csch, :sech, :coth)
    @eval @inline Base.$f(x::FiniteDual) = FiniteDual($f(Array(table(x))))
end

@inline ^(x::FiniteDual, y::Int) = FiniteDual(^(table(x), y))
@inline ^(x::FiniteDual, y::Number) = FiniteDual(^(Array(table(x)), y))
@inline ^(x::Number, y::FiniteDual) = FiniteDual(^(x, Array(table(y))))
@inline ^(x::FiniteDual, y::FiniteDual) = exp(y * log(x))

#################################
# Functions defined by branches #
#################################

struct CustomSign 
    result::FiniteDual
end

issign(::CustomSign) = true
issign(::FiniteDual) = false
issign(::Number) = false
issign(X::AbstractArray) = prod(map(x -> issign(x), X))

extract_sign(x::CustomSign) = x.result
extract_sign(x::FiniteDual) = x
extract_sign(x::Number) = x
extract_sign(X::AbstractArray) = map(x -> extract_sign(x), X)

"""
    custom_sign(x; fl::Function, fc::Function, fr::Function, a=0.0)

Customize the function definded with branches, defined by 
```math
F(x) = \\left\\{
\begin{aligned}
f_l(x), \\quad&{\rm if}\\,\\, x< a\\
f_c(x), \\quad&{\rm if}\\,\\, x= a\\
f_r(x), \\quad&{\rm if}\\,\\, x>a\\
\end{aligned}
\right.
```
"""

@inline function custom_sign(x::FiniteDual; fl::Function,
                             fc::Function, fr::Function, a=0.0)
    fd = sum(zip((<, ==, >), (fl, fc, fr))) do (pred, bfun)
        iszero(bfun(x)) ? zero(x) : bfun(x) * pred(x, a)
    end

    return CustomSign(fd)
end

function custom_sign(x::R; fl::FL, fc::FC, fr::FR, a=0.0) where {R<:Real,FL,FC,FR}
    y = x==a ? fc(x) : (x < a ? fl(x) : fr(x))
end

# standard sign function # 
#------------------------#

@inline Base.sign(x::FiniteDual) = custom_sign(x; fl=xl -> -1, fc=xc -> 0, fr=xr -> 1)

"""
    heaviside(x)

Return 1 if `x >= 0` and 0 otherwise.
"""

@inline heaviside(x::FiniteDual) = CustomSign(x >= 0)
@inline heaviside(x::Real) = x >= 0 ? 1 : 0
