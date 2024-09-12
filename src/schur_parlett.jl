#------------------------------------------------#
# Implemention of the schur parlett algorithm in #
# DOI: https://doi.org/10.1017/S0962492910000036 #
#------------------------------------------------#

# scalar function #
#-----------------#

function mat_fun_schur_parlett(f::Function, A::AbstractMatrix;
                               kwargs...)
    # Schur decomposition of A
    S = schur(A)
    if !istriu(S.T)
        S = Schur{Complex}(schur(A))
    end
    T, Z, Λ = S

    # For diagonalizable matrix A = Z*Λ*Z', 
    # f(A) = Z*f(Λ)*Z'
    isdiag(T) && return Z * diagm(fval(f, Λ)) * Z'

    # Split the eigenvalues into different blocks
    split = split_eigs_into_blocks(Λ; kwargs...)

    # If eigenvalues are distinct,
    # compute f(T) by standard Parlett recurrence
    split == 0 && return Z * parlett_recurence(f, T) * Z'

    # If there are close eigenvalues,
    # reorder the schur decomposition,
    # and compute f(T) by block Parlett recurrence
    reS, block_size = reorder_schur(S, split)
    return reS.Z * block_parlett_recurence(f, reS.T, block_size) * reS.Z'
end

# non-scalar function #
#---------------------#

function mat_fun_schur_parlett!(f!::Function, F::AbstractArray{V}, 
                                A::AbstractMatrix; 
                                kwargs...) where {V<:AbstractMatrix}
    
    # Schur decomposition of A
    S = schur(A)
    if !istriu(S.T)
        S = Schur{Complex}(schur(A))
    end
    T, Z, Λ = S

    # For diagonalizable matrix A = Z*Λ*Z', 
    # f(A) = Z*f(Λ)*Z'
    if isdiag(T)
        fval!(f!, F, FiniteDual(diagm(Λ)); kwargs...)
    else
        # Split the eigenvalues into different blocks
        split = split_eigs_into_blocks(Λ; kwargs...)

        # If eigenvalues are distinct,
        # compute f(T) by standard Parlett recurrence
        if split == 0
            parlett_recurence!(f!, F, T; kwargs...)
        else
            # If there are close eigenvalues,
            # reorder the schur decomposition,
            # and compute f(T) by block Parlett recurrence
            S, block_size = reorder_schur(S, split)
            T, Z, Λ = S
            block_parlett_recurence!(f!, F, T, block_size; kwargs...)
        end
    end

    map!(x -> Z * x * Z', F, F) 

    return F
end

#-------------------------------#
# Implemention of algorithm 4.2 #
#-------------------------------#

function parlett_recurence(f::Function, T::AbstractMatrix)
    @assert istriu(T)

    F = diagm(fval(f, diag(T)))
    for j = 2:size(F, 1), i = j-1:-1:1
        dt = T[i, i] - T[j, j]
        k = i+1:j-1
        F[i, j] = (T[i, j] * (F[i, i] - F[j, j]) / dt
                   +
                   (sum(F[i, k] .* T[k, j] .- T[i, k] .* F[k, j])) / dt)
    end

    return F
end

function parlett_recurence!(f!::Function, F::AbstractArray{V}, 
                            T::AbstractMatrix; 
                            kwargs...) where {V<:AbstractMatrix}
    @assert istriu(T)

    fval!(f!, F, FiniteDual(diag(T)); kwargs...)
    for j = 2:size(T, 1), i = j-1:-1:1
        dt = T[i, i] - T[j, j]
        k = i+1:j-1
        for (m, fm) in enumerate(F)
            F[m][i, j] = (T[i, j] * (fm[i, i] - fm[j, j]) / dt
                          +
                          (sum(fm[i, k] .* T[k, j] .- T[i, k] .* fm[k, j])) / dt)
        end
    end

    return F
end

#-------------------------------#
# Implemention of algorithm 4.3 #
#-------------------------------#

function block_parlett_recurence(f::Function, T::AbstractMatrix,
                                 block_size::Vector{<:Integer})
    @assert istriu(T)

    F = fill!(similar(T, typeof(f(T[1]))), 0)
    block = vcat(0, block_size)
    ind(i) = block[i]+1:block[i+1]
    for j = 1:length(block)-1
        indj = ind(j)
        lj = length(indj)
        Tjj = T[indj, indj]
        F[indj, indj] = mat_fun_atomic_block(f, FiniteDual(Tjj), Val(lj))
        for i = j-1:-1:1
            indi = ind(i)
            li = length(indi)
            Tii = T[indi, indi]
            F[indi, indi] = mat_fun_atomic_block(f, FiniteDual(Tii), Val(li))
            Y = F[indi, indi] * T[indi, indj] - T[indi, indj] * F[indj, indj]
            for k = i+1:j-1
                indk = ind(k)
                Y = Y + (F[indi, indk] * T[indk, indj]
                         -
                         T[indi, indk] * F[indk, indj])
            end

            # solve Tii*Fij + Fij*(-Tjj) + (-Y) = 0
            if li == 1 && lj == 1
                F[indi, indj] = Y ./ (Tii - Tjj)
            else
                F[indi, indj] = sylvester(Tii, -Tjj, -Y)
            end
        end
    end

    return F
end

function block_parlett_recurence!(f!::Function, F::AbstractArray{V},
                                  T::AbstractMatrix,
                                  block_size::Vector{Int64};
                                  kwargs...) where {V<:AbstractMatrix}

    block = vcat(0, block_size)
    ind(i) = block[i]+1:block[i+1]
    for j = 1:length(block)-1
        indj = ind(j)
        lj = length(indj)
        Tjj = T[indj, indj]
        Fjj = similar(F, typeof(Tjj))
        mat_fun_atomic_block!(f!, Fjj, FiniteDual(Tjj), Val(lj); kwargs...)
        for mj in eachindex(F)
            F[mj][indj, indj] = Fjj[mj]
        end
        for i = j-1:-1:1
            indi = ind(i)
            li = length(indi)
            Tii = T[indi, indi]
            Fii = similar(F, typeof(Tii))
            mat_fun_atomic_block!(f!, Fii, FiniteDual(Tii), Val(li); kwargs...)
            for mi in eachindex(F)
                F[mi][indi, indi] = Fii[mi]
                Y = F[mi][indi, indi] * T[indi, indj] - T[indi, indj] * F[mi][indj, indj]
                for k = i+1:j-1
                    indk = ind(k)
                    Y = Y + (F[mi][indi, indk] * T[indk, indj]
                             -
                             T[indi, indk] * F[mi][indk, indj])
                end

                # solve Tii*Fij + Fij*(-Tjj) + (-Y) = 0
                if li == 1 && lj == 1
                    F[mi][indi, indj] = Y ./ (Tii - Tjj)
                else
                    F[mi][indi, indj] = sylvester(Tii, -Tjj, -Y)
                end
            end
        end
    end

    return F
end

#------------------------------------------#
# Compute f(T),                            #
# for the eigenvalues of T are very close. #
# Directly use matrix functions in Julia   #
# instead of Taylor expansion.			   #
#------------------------------------------#

mat_fun_atomic_block(f, T::FiniteDual, ::Val{N}) where {N} = table(f(T))
mat_fun_atomic_block(f, T::AbstractMatrix, ::Val{1}) = f.(T)

function mat_fun_atomic_block!(f!::Function, F::AbstractArray, 
                               T::FiniteDual, ::Val{N}; 
                               kwargs...) where {N}
    fval!(f!, F, T; kwargs...)
end
function mat_fun_atomic_block!(f!::Function, F::AbstractArray{V}, 
                               T::AbstractMatrix, ::Val{1}; 
                               f!_T::Type=eltype(V), 
                               f!_S::Tuple=size(F),
                               kwargs...) where {V<:AbstractMatrix}
    tmp = Array{f!_T}(undef, f!_S)
    result = f!(tmp, T[1])
    @. F = ones(f!_T, 1, 1)
    result == nothing ? (@. F = tmp * F) : (@. F = result * F)
    return F
end

#------------------#
# Compute f.(eigs) #
#------------------#

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

function fval!(f!::Function, 
               F::AbstractArray{V},
               Λ::FiniteDual{T,N}; 
               f!_T::Type=eltype(V),
               f!_S::Tuple=size(F), 
               kwargs...) where {V<:AbstractMatrix,T,N}
    TT = FiniteDual{promote_type(T, f!_T),N}
    FD_tmp = Array{TT}(undef,f!_S)
    result = f!(FD_tmp, Λ)
    result == nothing ? map!(x -> table(x), F, FD_tmp) : 
                        map!(x -> table(x), F, result)
    return F
end

#------------------------------------------------------#
# Map λ to an integer q, s.t. λ ∈ S_q,                 #
# and satisfy the conditions: 						   #
# 1) min{|λ - μ|: λ ∈ S_p, μ ∈ S_q, p ≠ q} > δ.		   #
# 2) for S_p with |S_p| > 1, 						   #
#    ∀ λ ∈ S_p, ∃ μ ∈ S_p and μ ≠ λ, s.t. |λ - μ| ≤ δ. #
# For all S_q with |S_q| = 1, let q = 0.		       #
#------------------------------------------------------#

function split_eigs_into_blocks(eigs::Vector; δ::Real=0.1, kwargs...)
    # Sort the eigenvalues to quickly calculate distance
    sp = sortperm(eigs)
    eigs_sp = eigs[sp]
    N = length(eigs)
    dist = SVector{N - 1}(eigs_sp[2:end] - eigs_sp[1:end-1])

    minimum(dist) > δ && return 0

    split_map = Int[]
    split_pos = vcat(0, findall(x -> x > δ, dist), N) # find all split positions
    l = 1
    for i = 1:length(split_pos)-1
        dvec = dist[split_pos[i]+1:split_pos[i+1]-1]
        if length(dvec) == 0
            push!(split_map, 0)
        else
            append!(split_map, l * ones(Int, length(dvec) + 1))
            l += 1
        end
    end

    return split_map[sortperm(sp)]
end

#-------------------------------------------------------#
# Find the swap strategy that converts                  #
# an unordered sequence to a decreasing order sequence. #
# Note that in each step, 								#
# identical numbers are converted together.				#
# E.g., (1,0,2,1,2,3,0,3) -> (1,1,0,2,2,3,0,3) 			#
#       -> (2,2,1,1,0,3,0,3) -> (3,3,2,2,1,1,0,0).		#
#-------------------------------------------------------#

function swap_strategy(split)
    lmax = maximum(split)
    N = length(split)
    strategy = Vector{Bool}[]
    block_size = Int[]
    shift = 0
    for l = 1:lmax
        st_l = zeros(Bool, N)
        swap_index = findall(isequal(l), split)
        st_l[swap_index.+shift] .= 1
        push!(strategy, st_l)

        shift += length(swap_index)
        push!(block_size, shift)

        split = split[setdiff(1:length(split), swap_index)]
    end
    append!(block_size, collect(block_size[end]+1:N))

    return strategy, block_size
end

#---------------------------------# 
# Reorder the Schur decomposition #
#---------------------------------# 

function reorder_schur(S::Schur, split)
    strategy, block_size = swap_strategy(split)
    for st in strategy
        ordschur!(S, st)
    end
    return S, block_size
end
