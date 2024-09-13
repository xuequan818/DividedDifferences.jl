module SplitSwapTest

using Test
import DividedDifferences: split_eigs_into_blocks, swap_strategy

@testset "split points into blocks" begin
	# split distance
    dist = [0.05, 0.1]

    # distinct points
	distinct_points = collect(1:10)
	for δ in dist
        @test split_eigs_into_blocks(distinct_points) == 0
	end

	# close points
    close_points = [-0.2, -0.16, -0.1, -0.1, 0, 0.2, 0.36, 0.4]
    split = [[3, 3, 2, 2, 0, 0, 1, 1], [2, 2, 2, 2, 2, 0, 1, 1]]
    for (i, δ) in enumerate(dist)
        @test split_eigs_into_blocks(close_points; δ) == split[i]
    end
end

@testset "swap strategy" begin
    seq = [0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 3]
    order = [1, 3, 9, 7, 2, 4, 8, 6, 10, 5, 11]
    seq = seq[order]
	swap, bsize = swap_strategy(seq)	
	for sw in swap
		seq = vcat(seq[sw],seq[.!sw])
	end
    @test issorted(seq; rev=true)
    @test bsize == [0, 4, 6, 9, 10, 11]

    seq = [2, 2, 2, 1, 1, 0, 0]
    swap, bsize = swap_strategy(seq)
    @test swap == false
    @test bsize == [0, 3, 5, 6, 7]

    seq = [1, 1, 1]
    swap, bsize = swap_strategy(seq)
    @test swap == false
    @test bsize == [0, 3]
end

end