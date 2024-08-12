# DividedDifferences.jl

DividedDifferences computes **divided differences** by manipulating the **divided-difference tables** (see [![DOI](https://img.shields.io/badge/DOI-10.21105/jcon.00069-blue)](https://doi.org/10.1007/978-1-4020-6585-9_13)). While this method does not speed up execution compared to the naive approach (or even slower), its advantage is that the higher-order results are more accurate. In particular, when the points $x_0,x_1,\dots,x_n$ are very close or identical, the numerical error of the naive approach blows up, whereas this method ensures accurate convergence.

Since the divided-difference table for a function `f` is an upper triangular matrix, DividedDifferences supports the Julia matrix functions (or any object composed of them). Additionally, DividedDifferences also supports the special functions which are definded with branches, such as sign function `sign` and Heaviside step function `DividedDifferences.heaviside`. Users can also customize branching functions with the form
```math
F(x) = \left\{
\begin{aligned}
&f_l(x), \quad&{\rm if}\,\, x< a\\
&f_c(x), \quad&{\rm if}\,\, x= a\\
&f_r(x), \quad&{\rm if}\,\, x>a\\
\end{aligned}
\right.
```
by using `DividedDifferences.custom_sign(x; fl::Function, fc::Function, fr::Function, a)`.

Here are some simple examples showing how to use the package: 
```julia
julia> using DividedDifferences

julia> f(x) = DividedDifferences.custom_sign(x; fl=x->(x-1)/(x+1), fc=x->1.0, fr=x->0.0, a=1.); # returns a scalar

julia> x = [0.9, 0.9, 1.1]
3-element Vector{Float64}:
 0.9
 0.9
 1.1

julia> divided_difference(f, x) # returns the second order divided difference f[0.9, 0.9, 1.1]
-1.4542936288088644

julia> g(x) = [sin(x + 1) cos(x - 1); exp(x) x^3] # returns an array

julia> divided_difference(g, 1.0, 1.1) # returns the first order divided difference g[1.0, 1.1]
2×2 Matrix{Float64}:
 -0.460881  -0.0499583
  2.85884    3.31
```