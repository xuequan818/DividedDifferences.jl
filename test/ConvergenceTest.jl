module ConvergenceTest

using Test
using DividedDifferences
using DiffTests
using ForwardDiff

# compute the n-th order derivative
function derivative_nth_order(f::Function, x::R, n::Integer) where {R<:Real}
    D1 = x -> ForwardDiff.derivative(f, x)
    Ds = Function[D1]
    for i = 2:n
        push!(Ds, x -> ForwardDiff.derivative(Ds[i-1], x))
    end

    return Ds[end](x)
end

const x, N = 0.1, 4

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS
    for i = 2:N
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "$f" for f in DiffTests.NUMBER_TO_ARRAY_FUNCS
    for i = 2:N
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)

        out = similar(dd)
        out = divided_difference!(out, f, x * ones(i))
        @test isapprox(out, dd)
    end
end

@testset "$(f!)" for f! in DiffTests.INPLACE_NUMBER_TO_ARRAY_FUNCS
    m, n = 3, 2
    y = fill(0.0, m, n)
    f = x -> (tmp = similar(y, promote_type(eltype(y), typeof(x)), m, n); f!(tmp, x); tmp)

    for i = 2:N
        v = f(x)
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)

        fill!(y, 0.0)
        @test isapprox(divided_difference(f!, y, x * ones(i)), dd)
        @test isapprox(v, y)

        out = similar(dd)
        fill!(y, 0.0)
        divided_difference!(out, f!, y, x * ones(i))
        @test isapprox(out, dd)
        @test isapprox(v, y)
    end
end

@testset "spectial function defined by branches" begin
    f(x) = DividedDifferences.custom_sign(x; fl=xl -> exp(1 / (xl^2 + 1)), fc=xc -> 0, fr=xr -> cos(xr) - 1, a=3)
    for i = 2:N
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "heaviside step function" begin
    f(x) = DividedDifferences.heaviside(x)
    for i = 2:N
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "complex output" begin
    f(x) = (1 + im) * x
    for i = 2:N
        d = derivative_nth_order(f, x, i - 1)
        dd = divided_difference(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

end