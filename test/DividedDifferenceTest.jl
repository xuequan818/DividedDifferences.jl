module DividedDifferenceTest

using Test
using DividedDifferences
using DiffTests

# compute the divided difference by definition
function divided_difference_naive(f::F, x::Vector{R}) where {F,R<:Real}
    N = length(x)

    d1 = [(f(x[i]) - f(x[i+1])) / (x[i] - x[i+1]) for i = 1:N-1]
    dd = Vector[d1]
    for i = 2:N-1
        @views di = dd[i-1]
        push!(dd, [(di[j] - di[j+1]) / (x[j] - x[j+i]) for j = 1:N-i])
    end

    dd[end][1]
end

const X = [1, 2, 3, 4]

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS
    for i = 2:length(X)
		x = X[1:i]
		dd_naive = divided_difference_naive(f, x)
		dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "$f" for f in DiffTests.NUMBER_TO_ARRAY_FUNCS
    for i = 2:length(X)
		x = X[1:i]
		dd_naive = divided_difference_naive(f, x)
		dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)

        out = similar(dd_finitedual)
		out = divided_difference!(out, f, x)
        @test isapprox(out, dd_finitedual)
    end
end

@testset "$(f!)" for f! in DiffTests.INPLACE_NUMBER_TO_ARRAY_FUNCS
    m, n = 3, 2
    y = fill(0.0, m, n)
    f = x -> (tmp = similar(y, promote_type(eltype(y), typeof(x)), m, n); f!(tmp, x); tmp)

    for i = 2:length(X)
        x = X[1:i]
        val0 = f(x[1])
        dd_naive = divided_difference_naive(f, x)
        dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)

        fill!(y, 0.0)
        @test isapprox(divided_difference(f!, y, x), dd_finitedual)
        @test isapprox(val0, y)

        out = similar(dd_finitedual)
        fill!(y, 0.0)
        divided_difference!(out, f!, y, x)
        @test isapprox(out, dd_finitedual)
		@test isapprox(val0, y)
    end
end

@testset "spectial function defined by branches" begin
    f(x) = DividedDifferences.custom_sign(x - 3, xl -> exp(1 / (xl^2 + 1)), xc -> 0, xr -> cos(xr) - 1)
    for i = 2:length(X)
        x = X[1:i]
        dd_naive = divided_difference_naive(f, x)
        dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "heaviside step function" begin
	f(x) = DividedDifferences.heaviside(x-2)
	for i = 2:length(X)
        x = X[1:i]
        dd_naive = divided_difference_naive(f, x)
        dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "complex output" begin
    f(x) = (1+im)*x
    for i = 2:length(X)
        x = X[1:i]
        dd_naive = divided_difference_naive(f, x)
        dd_finitedual = divided_difference(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

end