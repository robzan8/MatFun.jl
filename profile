500x500 complex random matrix (separated eigenvalues, 1x1 blocks):
493 schur
89 schurparlett(T, Q)
	62 parlettrec
		34 trsyl
	24 return Q*F*Q'

500x500 complex random matrix ./ 60 (near eigenvalues, one big block):
494 schur
221 schurparlett(T, Q)
	218 atomicblock
		various matrix operations
		(TaylorSeries doesn't show up)

500x500 complex random matrix ./ 28 (22 blocks):
522 schur
246 schurparlett(T, Q)
	197 parlettrec
		172 atomicblock
