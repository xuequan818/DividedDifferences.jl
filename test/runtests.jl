using DividedDifferences
using Test

@testset "DividedDifferences.jl" begin
    t0 = time()
    @testset "Split and Swap" begin
        println("##### Testing split and swap strategy...")
        t = @elapsed include("SplitSwapTest.jl")
        println("##### done (took $t seconds).")
    end
    @testset "Schur Parlett Algorithm" begin
        println("##### Testing schur parlett algorithm functionality...")
        t = @elapsed include("DividedDifferenceTest.jl")
        println("##### done (took $t seconds).")
    end
    @testset "Convergence" begin
        println("##### Testing convergence of divided difference functionality...")
        t = @elapsed include("ConvergenceTest.jl")
        println("##### done (took $t seconds).")
    end
    println("##### Running all DividedDifferences tests took $(time() - t0) seconds.")
end
