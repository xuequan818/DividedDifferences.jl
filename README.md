[![Build Status](https://github.com/xuequan818/DividedDifferences.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/xuequan818/DividedDifferences.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/xuequan818/DividedDifferences.jl/graph/badge.svg?token=XNA1NK3RI8)](https://codecov.io/gh/xuequan818/DividedDifferences.jl)

# DividedDifferences.jl

DividedDifferences computes **divided differences** by manipulating the **divided-difference tables** (see [![DOI](https://img.shields.io/badge/DOI-10.21105/jcon.00069-blue)](https://doi.org/10.1007/978-1-4020-6585-9_13)). While this method does not speed up execution compared to the naive approach (or even slower), its advantage is that the higher-order results are more accurate. In particular, when the points $x_0,x_1,\dots,x_n$ are very close or identical, the numerical error of the naive approach blows up, whereas this method ensures accurate convergence.

Since the divided-difference table for a function `f` is an upper triangular matrix,  DividedDifferences supports the Julia matrix functions (or any object composed of them). Additionally, DividedDifferences also supports the special functions which are definded with branches, such as sign function `sign` and Heaviside step function `heaviside`. Users can also customize branch functions with the form
```math
F(x) = \left\{
\begin{aligned}
f_l(x), \quad&{\rm if}\,\, x< a\\
f_c(x), \quad&{\rm if}\,\, x= a\\
f_r(x), \quad&{\rm if}\,\, x>a\\
\end{aligned}
\right.
```
by using `custom_sign(x; fl::Function, fc::Function, fr::Function, a)`. Users need to set `div_diff(f, x; ill_test=false)` when compute the divided difference for branch functions. However, users should be aware that the results are not reliable when there are points on the discontinuities.

Here are some simple examples showing how to use the package: 
```julia
julia> using DividedDifferences

julia> f(x) = custom_sign(x; fl=x->(x-1)/(x+1), fc=x->1.0, fr=x->0.0, a=1.); # returns a scalar

julia> div_diff(f, 0.8, 0.8, 0.99, 1.01; ill_test=false) # returns the third order divided difference f[0.8, 0.8, 0.99, 1.01]
-5.486405741650227

julia> div_diff(heaviside, -0.1, 0.1, -0.1; ill_test=false) # returns the second order divided difference heaviside[-0.1, 0.1, -0.1]
25.0

julia> g(x) = [sin(x + 1) cos(x - 1); exp(x) x^3] # returns an array

julia> x = [1.0,1.1]
2-element Vector{Float64}:
 1.0
 1.1

julia> dd = div_diff(g, x) # returns the first order divided difference g[1.0, 1.1]
2×2 Matrix{Float64}:
 -0.460881  -0.0499583
  2.85884    3.31

julia> out = similar(dd);

julia> div_diff!(out, g, x); # store the divided difference in out

julia> out
2×2 Matrix{Float64}:
 -0.460881  -0.0499583
  2.85884    3.31
```

DividedDifferences can still deal with the non-scalar functions:
```julia
julia> function f!(y, x)                   # non-scalar function, store the result in y
	   fill!(y, zero(x))
           y[1] = x
           for i in 2:length(y)
               y[i] = exp(x) * sin(y[i-1])
           end
           return nothing
        end;

julia> y = zeros(3,2)
3×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
 0.0  0.0

julia> x = collect(0.9:0.1:1.1)
3-element Vector{Float64}:
 0.9
 1.0
 1.1

julia> dd = div_diff(f!, y, x) # returns the second order divided difference  
                                         # f![0.9, 1.0, 1.1] and stores f!(x[1]) in y
3×2 Matrix{Float64}:
   0.0       -3.58207
   1.46488  -30.2503
 -22.446    -89.8768

julia> y
3×2 Matrix{Float64}:
 0.9      1.82512
 1.92667  2.38049
 2.30549  1.69644

julia> out = similar(dd);

julia> fill!(y, 0.0);
3×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
 0.0  0.0

julia> div_diff!(out, f!, y, x); # stores divided difference in out and f!(x[1]) in y

julia> out
3×2 Matrix{Float64}:
   0.0       -3.58207
   1.46488  -30.2503
 -22.446    -89.8768

julia> y
3×2 Matrix{Float64}:
 0.9      1.82512
 1.92667  2.38049
 2.30549  1.69644
```