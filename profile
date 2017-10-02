1000x1000 real random matrix (separated eigenvalues):
time: 6s, samples: 3016.
2613 schur
341 recf
	192 trsyl

srand(666)
1000x1000 real random matrix / 19 (34 blocks):
time: 7.5s, samples: 3704.
2600 schur
352 reorder/trsen
717 recf
	343 taylorf
		174 gemm
	208 trsyl

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
