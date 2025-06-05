using StaticArrays
using LinearAlgebra
using IterTools

function mat_fun(f::Function, A::AbstractMatrix; 
				 splitting_criterion=0.1, tol=splitting_criterion)
    # Schur decomposition of A.
    S = schur(A)
    if !istriu(S.T)
        S = Schur{Complex}(schur(A))
    end
    T, Z, Λ = S

    # For diagonalizable matrix A = Z*Λ*Z', 
    # f(A) = Z*f(Λ)*Z'.
    isdiag(T) && return Z * diagm(fval(f, Λ)) * Z'

    # Split eigenvalues into blocks by `splitting_criterion`.
    split = points_splitting(Λ, splitting_criterion)

    # If all the eigenvalues are separated,
    # compute f(T) by standard Parlett recurrence.
    iszero(split) && return Z * parlett_recurrence(f, T) * Z'

    # If all the eigenvalues are close,
    # compute f(T) by `inner_mat_fun`.
    isone(split) && return Z * inner_mat_fun(f, T, tol) * Z'

    # If there are close eigenvalues,
    # reorder the schur decomposition,
    # and compute f(T) by block Parlett recurrence.
    reS, block_size = reorder_schur(S, split)
    reS.Z * block_parlett_recurrence(f, reS.T, block_size, tol) * reS.Z'
end

function parlett_recurrence(f::Function, T::AbstractMatrix)
    @assert istriu(T)

    F = diagm(fval(f, diag(T)))
    for j = 2:size(F, 1), i = j-1:-1:1
        dt = T[i, i] - T[j, j]
        k = i+1:j-1
        F[i, j] = (T[i, j] * (F[i, i] - F[j, j]) / dt + 
                   (sum(F[i, k] .* T[k, j] .- T[i, k] .* F[k, j])) / dt)
    end

    return F
end

function block_parlett_recurrence(f::Function, T::AbstractMatrix,
    							  block_size::Vector{<:Integer}, tol)
    @assert istriu(T)

    F = fill!(similar(T, typeof(f(T[1]))), 0)
    ind(i) = block_size[i]+1:block_size[i+1]
    for j = 1:length(block_size)-1
        indj = ind(j)
        Tjj = T[indj, indj]
        F[indj, indj] = inner_mat_fun(f, Tjj, tol)
        for i = j-1:-1:1
            indi = ind(i)
            Tii = T[indi, indi]
            F[indi, indi] = inner_mat_fun(f, Tii, tol)
            Y = F[indi, indi] * T[indi, indj] - T[indi, indj] * F[indj, indj]
            for k = i+1:j-1
                indk = ind(k)
                Y = Y + (F[indi, indk] * T[indk, indj] - 
                         T[indi, indk] * F[indk, indj])
            end

            # solve Tii*Fij + Fij*(-Tjj) + (-Y) = 0
            if length(indi) == 1 && length(indj) == 1
                F[indi, indj] = Y ./ (Tii - Tjj)
            else
                F[indi, indj] = sylvester(Tii, -Tjj, -Y)
            end
        end
    end

    return F
end

# Compute f.(eigs) 
function fval(f::Function, Λ::Vector)
    try
        f.(Λ)
    catch e
        if isa(e, DomainError)
            f.(complex.(Λ))
        else
            throw(e)
        end
    end
end

# Inner matrix function thast computes f(T)
function inner_mat_fun(f::Function, A::AbstractMatrix, tol::Real)
    N = size(A,1)
    N == 1 && return f.(A)
					   
	if native(f) # check if f is a Julia native matrix function
		return f(A)
	elseif spread(A) > tol # gurantee the spread of block is smaller than `tol`
		return mat_fun(f, A; splitting_criterion=max(tol/N,cbrt(eps())), tol)
	else
        return taylor_mat_fun(f, A; max_iter=min(10N,100))
	end
end

const NATIVE_TEST_MAT = SDiagonal{2}(I)

function native(f::Function)
    try
        f(NATIVE_TEST_MAT)
        true
    catch
        false
    end
end

function spread(A::AbstractMatrix)
    dA = diag(A)
    maximum(dA) - minimum(dA)
end
