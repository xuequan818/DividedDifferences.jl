module SchurParlettTest

using Test
using DiffTests
using DividedDifferences
import DividedDifferences: table, mat_fun_schur_parlett, mat_fun_schur_parlett!

const X1, X2 = [3, 4, 5], [1.1, 1.3, 1.1, 1.4]

FD = map(x->FiniteDual(x), [X1, X2])

@testset "$f" for f in DiffTests.NUMBER_TO_NUMBER_FUNCS[1:3]
    for fd in FD
		@test isapprox(table(f(fd)), mat_fun_schur_parlett(f, table(fd)))
	end  
end

@testset "$f" for f in DiffTests.NUMBER_TO_ARRAY_FUNCS
    for fd in FD
        mat_op = Array.(table.(f(fd)))
        mat_sp = zero.(mat_op)
        f!(y, x) = copy!(y, f(x))
        mat_fun_schur_parlett!(f!, mat_sp, table(fd))
        @test isapprox(mat_op, mat_sp)
    end
end

@testset "$(f!)" for f! in DiffTests.INPLACE_NUMBER_TO_ARRAY_FUNCS
    m, n = 3, 2
    y = fill(0.0, m, n)
    g = x -> (tmp = similar(y, promote_type(eltype(y), typeof(x)), m, n); f!(tmp, x); tmp)

    for fd in FD
        mat_op = Array.(table.(g(fd)))
        mat_sp = zero.(mat_op)
        g!(y, x) = copy!(y, g(x))
        mat_fun_schur_parlett!(g!, mat_sp, table(fd))
        @test isapprox(mat_op, mat_sp)
    end

    for fd in FD
        yfd = similar(y, typeof(fd))
        f!(yfd, fd)
        mat_op = Array.(table.(yfd))
        mat_sp = zero.(mat_op)
        mat_fun_schur_parlett!(f!, mat_sp, table(fd))
        @test isapprox(mat_op, mat_sp)
    end
end

@testset "complex output" begin
    f(x) = (1 + im) * x
    for fd in FD
        @test isapprox(table(f(fd)), mat_fun_schur_parlett(f, table(fd)))
    end
end

end