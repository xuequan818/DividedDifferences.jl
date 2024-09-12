module DividedDifferenceTest

using Test
using DividedDifferences
using DiffTests

# compute the divided difference by definition
function div_diff_naive(f::F, x::Vector{R}) where {F,R<:Real}
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
		dd_naive = div_diff_naive(f, x)
		dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "$f" for f in DiffTests.NUMBER_TO_ARRAY_FUNCS
    for i = 2:length(X)
		x = X[1:i]
		dd_naive = div_diff_naive(f, x)
		dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)

        out = similar(dd_finitedual)
		out = div_diff!(out, f, x)
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
        dd_naive = div_diff_naive(f, x)
        dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)

        fill!(y, 0.0)
        @test isapprox(div_diff(f!, y, x), dd_finitedual)
        @test isapprox(val0, y)

        out = similar(dd_finitedual)
        fill!(y, 0.0)
        div_diff!(out, f!, y, x)
        @test isapprox(out, dd_finitedual)
		@test isapprox(val0, y)
    end
end

@testset "spectial function defined by branches" begin
    f(x) = custom_sign(x; fl=xl -> exp(1 / (xl^2 + 1)), fc=xc -> 0, fr=xr -> cos(xr) - 1, a=2.5)
    for i = 2:length(X)
        x = X[1:i]
        dd_naive = div_diff_naive(f, x)
        dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "heaviside step function" begin
	f(x) = heaviside(x-2.5)
	for i = 2:length(X)
        x = X[1:i]
        dd_naive = div_diff_naive(f, x)
        dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

@testset "complex output" begin
    f(x) = (1+im)*x
    for i = 2:length(X)
        x = X[1:i]
        dd_naive = div_diff_naive(f, x)
        dd_finitedual = div_diff(f, x)
        @test isapprox(dd_naive, dd_finitedual)
    end
end

end