module ConvergenceTest

using Test
using DividedDifferences
using DiffTests
using ForwardDiff

# compute the n-th order derivative
function nth_derivative(f::Function, x::Number, n::Integer)
    iszero(n) ? f(x) : ForwardDiff.derivative(x -> nth_derivative(f, x, n-1), x)
end

const x, N = 0.1, 4

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS
    for i = 2:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "$f" for f in DiffTests.NUMBER_TO_ARRAY_FUNCS
    for i = 2:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)

        out = similar(dd)
        out = div_diff!(out, f, x * ones(i))
        @test isapprox(out, dd)
    end
end

@testset "$(f!)" for f! in DiffTests.INPLACE_NUMBER_TO_ARRAY_FUNCS
    m, n = 3, 2
    y = fill(0.0, m, n)
    f = x -> (tmp = similar(y, promote_type(eltype(y), typeof(x)), m, n); f!(tmp, x); tmp)

    for i = 2:N
        v = f(x)
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)

        fill!(y, 0.0)
        @test isapprox(div_diff(f!, y, x * ones(i)), dd)
        @test isapprox(v, y)

        out = similar(dd)
        fill!(y, 0.0)
        div_diff!(out, f!, y, x * ones(i))
        @test isapprox(out, dd)
        @test isapprox(v, y)
    end
end

@testset "spectial function defined by branches" begin
    f(x) = custom_sign(x; fl=xl -> exp(1 / (xl^2 + 1)), fc=xc -> 0, fr=xr -> cos(xr) - 1, a=3)
    for i = 2:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "heaviside step function" begin
    f(x) = heaviside(x)
    for i = 2:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

@testset "complex output" begin
    f(x) = (1 + im) * x
    for i = 2:N
        d = nth_derivative(f, x, i - 1)
        dd = div_diff(f, x * ones(i))
        @test isapprox(d / prod(1:i-1), dd)
    end
end

end