using IterTools
#-----------------------------#
# Reorder Schur decomposition #
#-----------------------------#

# Map λ to an integer q, s.t. λ ∈ S_q,                 
# and satisfy the conditions: 						   
# 1) min{|λ - μ|: λ ∈ S_p, μ ∈ S_q, p ≠ q} > δ.		   
# 2) for S_p with |S_p| > 1, 						   
#    ∀ λ ∈ S_p, ∃ μ ∈ S_p and μ ≠ λ, s.t. |λ - μ| ≤ δ. 
# For all S_q with |S_q| = 1, let q = 0.	

# Real points
function points_splitting(pts::Vector{<:Real}, δ::Real)
    @assert length(pts) > 1

    # Sort the eigenvalues to quickly calculate distance
    sp = sortperm(pts; rev=true) # real pts
    pts_sp = pts[sp]
    N = length(pts)
    dist = pts_sp[1:end-1] - pts_sp[2:end]

    # All points are single point sets.
    minimum(dist) > δ && return 0

    # All points are in one set.
    maximum(dist) <= δ && return 1

    split_map = Int[]
    split_pos = vcat(0, findall(x -> x > δ, dist), N) # find all split positions
    dvec = split_pos[2:end] - split_pos[1:end-1]
    l = 1
    for di in dvec
        if di == 1
            push!(split_map, 0)
        else
            append!(split_map, repeat([l],di))
            l += 1
        end
    end

    return split_map[sortperm(sp)]
end

# Complex points
function points_splitting(pts::Vector{<:ComplexF64}, δ::Real)
    #dmat = map(x -> norm(x[1] - x[2]), Iterators.product(pts, pts))
    #StrictUpperVec(dmat)

    @time dist, index = sep_vec(pts)

    # All points are single point sets.
    minimum(dist) > δ && return 0

    # All points are in one set.
    maximum(dist) <= δ && return 1

    sp = sortperm(dist)
    split_pos = findlast(x -> x <= δ, sort!(dist))
    CLOSE_PTS = index[sp,:][1:split_pos,:]

    MAP = Vector{Int}[]
    @views for k = 1:split_pos
        Sk = CLOSE_PTS[k,:]
		indk = Int[]
		for (l, Sl) in enumerate(MAP)
			if Sk[1] in Sl || Sk[2] in Sl
				push!(indk, l)
			end
		end
      
		if length(indk) == 0
			push!(MAP, Sk)
		else
  			append!(MAP[indk[1]], Sk)
     		if length(indk) > 1            
				for i in indk[2:end]
					append!(MAP[indk[1]], MAP[i])
				end
				deleteat!(MAP,indk[2:end])
			end
		end
	end
	unique!.(MAP)

    # All points are in one set.
    length(MAP[1]) == length(pts) && return 1
                 
    split_map = zeros(Int, length(pts))
    ℓ = 1
    for S in MAP
        split_map[S] .= ℓ
        ℓ += 1
    end

    return split_map
end

function sep_vec(pts::Vector{T}) where {T<:Complex}
    N = length(pts)
    lvec = cld(N^2 - N, 2)
    dvec = Vector{T}(undef, lvec)
    ind = zeros(Int, lvec, 2)

    l = 0
    @views for i = 1:N, j = i+1:N
        l += 1
        dvec[l] = pts[i] - pts[j]
        ind[l,1] = i
        ind[l,2] = j
    end
    @assert l == lvec

    norm.(dvec), ind
end
function StrictUpperVec(A::AbstractMatrix)
    n = size(A, 1)
    Avec = Vector{eltype(A)}(undef, Int(n * (n - 1) / 2))
    ind = zeros(Int, length(Avec), 2)
    l = 1
    for i = 1:n, j = i+1:n
        Avec[l] = A[i,j]
        ind[l,1] = i
        ind[l,2] = j
        l += 1
    end
    Avec, ind
end

# Find the swap strategy that converts                  
# an unordered sequence to a decreasing order sequence. 
# Note that in each step, 								
# identical numbers are converted together.				
# E.g., (1,0,2,1,2,3,0,3) -> (1,1,0,2,2,3,0,3) 			
#       -> (2,2,1,1,0,3,0,3) -> (3,3,2,2,1,1,0,0).		
function swap_strategy(split)
    lmax = maximum(split)
    N = length(split)

    if issorted(split; rev=true)
        count = findlast(isequal(lmax), split)
        block_size = [0, count]
        for l = lmax:-1:2
            count += findlast(isequal(l - 1), split) - findlast(isequal(l), split)
            push!(block_size, count)
        end
        append!(block_size, collect(block_size[end]+1:N))
        return false, block_size
    end

    strategy = Vector{Bool}[]
    count = N - length(findall(iszero, split))
    block_size = collect(count:N)
    shift = 0
    for l = 1:lmax
        st_l = zeros(Bool, N)
        swap_index = findall(isequal(l), split)
        st_l[swap_index.+shift] .= 1
        push!(strategy, st_l)

        sl = length(swap_index)
        shift += sl
        count -= sl
        pushfirst!(block_size, count)

        split = split[setdiff(1:length(split), swap_index)]
    end

    return strategy, block_size
end

function reorder_schur(S::Schur, split)
    strategy, block_size = swap_strategy(split)
    if strategy != false
        for st in strategy
            ordschur!(S, st)
        end
    end
    return S, block_size
end
