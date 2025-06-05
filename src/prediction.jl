# View the prediction as piecewise constant function, 
# where the true corresponding to 1, and false corresponding to 0.
# The prediction function `f` acts on a matrix `A` could be computed by 
# solving the Sylvester equation for 2×2 block matrices. 
# Compute `f(A)` by the Schur decomposition `f(A) = Z*f(T)*Z'`.
# First, reorder and divide `T` into a 2×2 matrix of upper triangular blocks, 
# where one of the diagonal blocks has all eigenvalues 
# meeting the prediction and the other does not.
# Utilize f(A) = \sum_{k=0}^\infty f^{(k)}(\sigma)/k! (A-\sigma I), 
# the corresponding diagonal parts of `F := f(T)` are 
# the unit and zero matrices. By the commutativity relation FT = TF, 
# there is the Sylvester equation: `T11*F12 - F12*T22 + (T12*F22-F11*T12) = 0`.

const BINARY_PREDICATES = Symbol[:isequal, :isless, :<, :>, :(==), :(!=), :(<=), :(>=)]

function pred_mat_fun(pred_val::BitVector, S::Schur)
	N = size(S.T,1)
	TT = eltype(S)

    sum(pred_val) == 0 && return zero(S.T)
    sum(pred_val) == N && return Matrix{TT}(I, N, N)

	# reorder the schur
    reS, block_size = reorder_schur(S, Int.(pred_val) .+ 1)
    T, Z, _ = reS

    @assert length(block_size) == 3

	fx = TT.(sort(pred_val; rev=true))
    F = diagm(fx)
    i, j = [block_size[l]+1:block_size[l+1] for l = 1:2]
    Y = T[i, j] * F[j, j] - F[i, i] * T[i, j]

    # solve Tii*Fij - Fij*Tjj + Y = 0
    if length(i) > 1 || length(j) > 1
        F[i, j] = sylvester(T[i, i], -T[j, j], Y)
    else
        F[i, j] = Y ./ (T[i, i] - T[j, j])
    end

    return Z * F * Z'
end

for pred in BINARY_PREDICATES
    @eval @inline function mat_fun($pred, A::AbstractMatrix, y::Number)
		S = schur(A)
		pred_val = $pred.(S.values, y)
		return pred_mat_fun(pred_val, S)
	end

    @eval @inline function mat_fun($pred, x::Number, A::AbstractMatrix)
        S = schur(A)
        pred_val = $pred.(x, S.values)
        return pred_mat_fun(pred_val, S)
    end
end

#function mat_fun(f::typeof(sign), A::AbstractMatrix; kwargs...)
#	mat_fun(!=, A, zero(eltype(A))) * mat_fun(x->x/abs(x), A; kwargs...)
#end