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
`ratkrylovf(f, A, b, mmax=100, tol=1e-13)`<br />
The poles for the rational Krylov decomposition are given by the AAA rational approximation of f, with parameters mmax and tol. mmax effectively bounds the size of the Krylov space. If you want to specify the poles manually, do:<br />
`ratkrylovf(f, A, b, p)`<br />
The rational Krylov decomposition can be computed with<br />
`V, K, H = ratkrylov(A, b, p)`<br />
Again, when A is real, `f(conj(x)) == conj(f(x))` is assumed.
