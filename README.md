## MatFun
This package provides methods for computing matrix functions. It currently works with Float64 precision only.

## Schur-Parlett
Schur-Parlett can be used to compute f(A) for dense matrices. It uses higher-order automatic differentiation (TaylorSeries.jl, specifically). It can be called with:
`schurparlett(f, A)`
If you want to reuse the Schur decomposition of A, you can also use:
`schurparlett(f, T, Q, vals)`
The algorithm has a couple of performance improvements, compared to the one described in the paper. The Parlett recurrence is implemented in a cache-oblivious fashion and the algorithm works mostly in real arithmetic, when A is real (in which case, `f(conj(x)) == conj(f(x))` is assumed).
