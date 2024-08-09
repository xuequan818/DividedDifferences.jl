module DividedDifferences

using StaticArrays
using LinearAlgebra
using Base: require_one_based_indexing

export FiniteDual
export divided_difference
export divided_difference!
include("finitedual.jl")
include("divided_difference.jl")

end # module DividedDifferences
