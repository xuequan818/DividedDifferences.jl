module DividedDifferences

using StaticArrays
using LinearAlgebra
using Base: require_one_based_indexing

export FiniteDual
export custom_sign
export heaviside
export div_diff
export div_diff!
include("finitedual.jl")
include("schur_parlett.jl")
include("divided_difference.jl")

end # module DividedDifferences
