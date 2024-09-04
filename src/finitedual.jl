import Base: +, -, *, /, ^, >, <, ==, >=, <=

##############
# FiniteDual #
##############

struct FiniteDual{T,N} <: Number
    DD_table::SMatrix{N, N, T}
end

################
# Constructors #
################

@inline function FiniteDual(x::AbstractMatrix{T}) where {T}
	if !istriu(x)
		error("Cannot create a divided difference table")
    end
	
    N = size(x,1)
    FiniteDual(SMatrix{N,N,T}(x))
end

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
@inline Base.length(::FiniteDual{T,N}) where {T,N} = N

@inline Base.values(x::FiniteDual) = diag(x.DD_table)
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
    if isconcretetype(R) # issue #322
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

@inline Base.convert(::Type{FiniteDual{T,N}}, x::FiniteDual{V,N}) where {T,V,N} = FiniteDual(T.(x.DD_table))
@inline Base.convert(::Type{FiniteDual{T,N}}, x::Number) where {T,N} = FiniteDual{N}(T(x))
@inline Base.convert(::Type{FD}, x::FD) where {FD<:FiniteDual} = x

#####################
# Generic Functions #
#####################

@inline Base.zero(x::FiniteDual) = zero(typeof(x))
@inline Base.zero(::Type{FiniteDual{T,N}}) where {T,N} = FiniteDual{N}(T(0))

@inline Base.one(x::FiniteDual) = one(typeof(x))
@inline Base.one(::Type{FiniteDual{T,N}}) where {T,N} = FiniteDual{N}(T(1))

# Predicates #
#------------#

isconstant(x::FiniteDual) = values(x) == value(x,1) * I
Base.iszero(x::FiniteDual) = iszero(x.DD_table)

for pred in [:>, :<, :(==), :>=, :<=]
    @eval @inline $pred(x::FiniteDual, y::Number) = 
                   build_pred_table(x, $pred.(diag(x.DD_table), y))
    @eval @inline $pred(x::Number, y::FiniteDual) = 
                   build_pred_table(y, $pred.(x, diag(y.DD_table)))
end

@inline function build_pred_table(x::FiniteDual{T,N},
                                  pred_val::SVector{N,Bool}) where {T,N}
    fx_val = values(x)

    if issorted(fx_val)
        return build_pred_table_in_sort(fx_val, pred_val)
    end

    pred_table = MMatrix{N,N}(diagm(float(T).(pred_val)))
    for k = 1:N-1
        indk = diagind(pred_table, k)
        for (l, indl) in enumerate(indk)
            pts = fx_val[l:l+k]
            sp = sortperm(pts)
            table_in_sort = build_pred_table_in_sort(pts[sp], pred_val[l:l+k][sp])
            pred_table[indl] = extract_DD(table_in_sort)
        end
    end

    return FiniteDual(pred_table)
end

function build_pred_table_in_sort(fx_val::AbstractVector,
                                  pred_val::AbstractVector)
    @assert issorted(fx_val)

    T = eltype(fx_val)
    N = length(fx_val)

    if sum(pred_val) in (0, N)
        return FiniteDual{N}(T(pred_val[1]))
    else
        pred_table = MMatrix{N,N}(diagm(float(T).(pred_val)))
        for k = 1:N-1
            indk = diagind(pred_table, k)
            valk = diag(pred_table, k - 1)
            valkdiff = valk[2:end] - valk[1:end-1]
            ptsdiff = fx_val[k+1:N] - fx_val[1:N-k]
            diffBit = (valkdiff .!= 0)
            pred_table[indk[diffBit]] = valkdiff[diffBit] ./ ptsdiff[diffBit]
        end
        return FiniteDual(pred_table)
    end
end

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
		FiniteDual($op(x.DD_table,y.DD_table))
	end
    @eval @inline $op(x::FiniteDual{Tx,Nx}, y::Number) where {Tx,Nx} = FiniteDual($op(x.DD_table, SDiagonal{Nx}(y * I)))
    @eval @inline $op(x::Number, y::FiniteDual{Ty,Ny}) where {Ty,Ny} = FiniteDual($op(SDiagonal{Ny}(x * I), y.DD_table))
end
@inline -(x::FiniteDual) = FiniteDual(-x.DD_table)

# (*)/(/) #
#---------#

for op in (:*, :/)
    @eval @inline function $op(x::FiniteDual, y::FiniteDual)
        @assert length(x) == length(y)
        FiniteDual($op(x.DD_table, y.DD_table))
    end
    @eval @inline $op(x::FiniteDual, y::Number) = FiniteDual($op(x.DD_table, y))
end
@inline *(x::Number, y::FiniteDual) = FiniteDual(x * y.DD_table)
@inline /(x::Number, y::FiniteDual{T,N}) where {T,N} = FiniteDual(\(y.DD_table, SDiagonal{N}(T(x) * I)))

# matrix functions #
#------------------#

for f in (:exp, :cis, :sqrt,
          :asin, :acos, :atan,
          :acsc, :asec, :acot,
          :asinh, :acosh, :atanh,
          :acsch, :asech, :acoth)
    @eval @inline Base.$f(x::FiniteDual) = FiniteDual($f(x.DD_table))
end

for f in (:log, :sin, :cos, :tan,
          :csc, :sec, :cot,
          :sinh, :cosh, :tanh,
          :csch, :sech, :coth)
    @eval @inline Base.$f(x::FiniteDual) = FiniteDual($f(Array(x.DD_table)))
end

@inline ^(x::FiniteDual, y::Int) = FiniteDual(^(x.DD_table, y))
@inline ^(x::FiniteDual, y::Number) = FiniteDual(^(Array(x.DD_table), y))
@inline ^(x::Number, y::FiniteDual) = FiniteDual(^(x, Array(y.DD_table)))
@inline ^(x::FiniteDual, y::FiniteDual) = exp(y * log(x))

# custom functions defined by sign rule #
# custom_sign(x) = fl(x) if x < a;      #
#                  fc(x) if x = a;      #
#                  fr(x) if x > a.      #
#---------------------------------------#

@inline function custom_sign(x::FiniteDual{T}; 
                             fl::FL, fc::FC, fr::FR,
                             a=0.0) where {T,FL,FC,FR}
    fd = sum(zip((<, ==, >), (fl, fc, fr))) do (pred, bfun)
        iszero(bfun(x)) ? zero(x) : bfun(x) * pred(x, T(a))
    end

    return fd
end
function custom_sign(x::R; fl::FL, fc::FC, fr::FR, a=0.0) where {R<:Real,FL,FC,FR}
    y = x==a ? fc(x) : (x < a ? fl(x) : fr(x))
end

# standard sign function # 
#------------------------#

@inline Base.sign(x::FiniteDual) = custom_sign(x; fl=xl -> -1, fc=xc -> 0, fr=xr -> 1)

# heaviside step function #
#-------------------------#

@inline heaviside(x) = custom_sign(x; fl=xl -> 0, fc=xc -> 1, fr=xr -> 1)

# Accurately compute 1/(1+exp(x)).                            #
# Rewritten by log-sum-exp trick to avoid overflow/underflow. #
# First write 1 / (1+exp(x)) = exp(-log(1+exp(x))),           #
# then log(1+exp(x)) = x + log(1+exp(x)) if x ≥ 0,            #
#                    = log(1+exp(x))     if x < 0.            #
# When act on FiniteDual and                                  #
# the `x` ponits are not on the same branch,                  #
# use the average of two branches:                            #
# 0.5 * x + log(exp(-0.5 * x) + exp(0.5 * x)).                   #
# In this case, it is valid for `abs(x[i]) ≤ 1000`.           #
#-------------------------------------------------------------#

@inline function invexp1p(x::FiniteDual) 
    fx_bit = values(x) .>= 0
    if sum(fx_bit) == length(x)
        return exp(-(x + log(1+exp(-x))))
    elseif sum(fx_bit) == 0
        return exp(-log(1+exp(x)))
    else
        return exp(-(0.5 * x + log(exp(-0.5 * x) + exp(0.5 * x))))
    end
end

@inline invexp1p(x::Real) = custom_sign(x; fl=xl -> exp(-log(1 + exp(xl))), fc=xc -> 0.5, fr=xr -> exp(-(xr + log(1 + exp(-xr)))))
    