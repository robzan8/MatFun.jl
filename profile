500x500 complex random matrix (separated eigenvalues, 1x1 blocks):
489 schur
13287 schurparlett(T, Q)
	12991 C += F[I,K]*T[K,J] - T[I,K]*F[K,J]
		5159 abstractarray getindex
		4900 matmul

500x500 complex random matrix ./ 60 (near eigenvalues, one big block):
677 schur
259 schurparlett(T, Q)
	246 atomicblock
		various matrix operations
		(TaylorSeries doesn't show up)

500x500 complex random matrix ./ 28 (22 blocks):
505 schur
429 schurparlett(T, Q)
	252 atomicblock (various matrix operations)
	31 C = F[I,I]*T[I,J] - T[I,J]*F[J,J] (getindex)
	68 trsyl
