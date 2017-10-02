1000x1000 complex random matrix (separated eigenvalues, 1x1 blocks):
time: 16s, samples: 6924.
4328 schur
2596 schurparlett(T, Q)
	1921 C += F[I,K]*T[K,J] - T[I,K]*F[K,J]
		1614 matmul
	475 trsyl getindex
	197 trsyl

srand(666)
1000x1000 complex random matrix ./ 19 (60 blocks):
time: 16s, samples: 7213.
4413 schur
2800 schurparlett(T, Q)
	1257 trexc
	451 atomicblock (various matrix operations)
	384 C = F[I,I]*T[I,J] - T[I,J]*F[J,J]
	676 trsyl getindex
	584 trsyl

500x500 complex random matrix ./ 60 (near eigenvalues, one big block):
677 schur
259 schurparlett(T, Q)
	246 atomicblock
		various matrix operations
		(TaylorSeries doesn't show up)
