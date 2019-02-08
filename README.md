(Make sure to take a look at the paper! https://github.com/pinkgopher/MatFun.jl/blob/master/docs/thesis.pdf )

## MatFun
This package provides methods for computing matrix functions. It currently works with Float64 precision only.

## Schur-Parlett
Schur-Parlett can be used to compute f(A) for dense matrices. It uses higher-order automatic differentiation (TaylorSeries.jl, specifically). It can be called with:<br />
`schurparlett(f, A)`<br />
If you want to reuse the Schur decomposition of A, you can also use:<br />
`schurparlett(f, T, Q, vals)`<br />
The algorithm has a couple of performance improvements, compared to the one described in the paper. The Parlett recurrence is implemented in a cache-oblivious fashion and the algorithm works mostly in real arithmetic, when A is real (in which case, `f(conj(x)) == conj(f(x))` is assumed).

## Rational Krylov
Rational Krylov can be used to compute f(A)*b for sparse matrices:<br />
`ratkrylov(f, A, b, mmax=100, tol=1e-13, Z=Vector{Complex128}(0))`<br />
The poles for the rational Krylov decomposition are given by the AAA rational approximation of f, with parameters mmax, tol and Z. mmax effectively bounds the size of the Krylov space. When the sample set Z is not provided, f will be sampled on the 0-centered disk with radius min(norm(A, 1), norm(A, Inf), vecnorm(A)). If you want to specify the poles manually, do:<br />
`ratkrylov(f, A, b, p)`<br />
The rational Krylov decomposition can be computed with:<br />
`V, K, H = ratkrylov(A, b, p)`<br />
Again, when A is real, operations are done in real arithmetic and `f(conj(x)) == conj(f(x))` is assumed.

## AAA Rational Approximation
The package also includes the AAA algorithm for rational approximation:<br />
`r, pol, res, zer, z, f, w, errvec = aaa(func, Z, tol=1e-13, mmax=100)`<br />
The algorithm is as described in the original paper.

## Installation
Do:<br />
`Pkg.add("TaylorSeries")`<br />
`Pkg.clone("https://github.com/pinkgopher/MatFun.jl.git")`<br />
Before using the package, make sure to test it with:<br />
`Pkg.test("MatFun")`

## Enjoy!
