#=
This file implements the algorithm described in:
"A Schur-Parlett Algorithm for Computing Matrix Functions"
(Davies and Higham, 2003)
=#

#=
blockpattern groups eigenvalues in sets, implementing Algorithm 4.1 of the paper.
vals is the vector of eigenvalues, δ is the tolerance.
The function returns S and p. S is the pattern, where S[i] = s means that
vals[i] has been assigned to set s. p is the number of sets identified.
=#
function blockpattern(vals::Vector{C}, δ::Float64) where {C}
	S = fill(-1, length(vals)) # -1 means unassigned
	p = 0
	for i = 1:length(vals)
		λ = vals[i]
		if S[i] == -1
			p += 1
			# assign λ to set p:
			for k = i:length(vals)
				if vals[k] == λ
					S[k] = p
				end
			end
		end
		Sλ = S[i]
		for j = i+1:length(vals)
			μ = vals[j]
			Sμ = S[j]
			if Sμ != Sλ && abs(λ-μ) <= δ
				if Sμ == -1
					# assign μ to set Sλ:
					for k = j:length(vals)
						if vals[k] == μ
							S[k] = Sλ
						end
					end
				else
					# merge sets Sλ and Sμ:
					Smin, Smax = minmax(Sλ, Sμ)
					for k = 1:length(vals)
						if S[k] == Smax
							S[k] = Smin
						elseif S[k] > Smax
							S[k] -= 1
						end
					end
					Sλ = Smin
					p -= 1
				end
			end
		end
	end
	return S, p
end

#=
reorder reorders the complex Schur decomposition (T, Q) according to pattern S,
using the swapping strategy described in Algorithm 4.2.
The entries of S are also reordered together with the corresponding eigenvalues.
The function returns vector blocksize: blocksize[i] is the size of the i-th
leading block on T's diagonal (after the reordering).
=#
function reorder!(T::Matrix{C}, Q::Matrix{C}, S::Vector{Int64}, p::Int64) where {C<:Complex}
	blocksize = zeros(Int64, p)
	# for each set, calculate its mean position in S:
	pos = zeros(Float64, p)
	count = zeros(Int64, p)
	for i = 1:length(S)
		set = S[i]
		pos[set] += i
		count[set] += 1
	end
	pos ./= count

	ilst = 1
	for set = 1:p
		minset = indmin(pos)
		for ifst = ilst:length(S)
			if S[ifst] == minset
				if ifst != ilst
					LAPACK.trexc!('V', ifst, ilst, T, Q)
					S[ilst], S[ilst+1:ifst] = S[ifst], S[ilst:ifst-1]
				end
				ilst += 1
				blocksize[set] += 1
			end
		end
		pos[minset] = Inf
	end
	return blocksize
end

#=
atomicblock computes f(T) using Taylor as described in Algorithm 2.6.
=#
function atomicblock(f::Func, T::Matrix{C}) where {Func, C}
	n = size(T, 1)
	if n == 1
		return f.(T)
	end

	maxiter = 300
	lookahead = 10
	σ = mean(diag(T))
	tay = f(σ + Taylor1(typeof(σ), 20+lookahead))
	M = UpperTriangular(T - σ*eye(T))
	P = M
	F = UpperTriangular(tay.coeffs[1]*eye(T))
	for k = 1:maxiter
		needorder = k + lookahead
		@assert needorder <= tay.order + 1
		if needorder > tay.order
			tay = f(σ + Taylor1(typeof(σ), min(maxiter+lookahead, tay.order*2)))
		end

		Term = tay.coeffs[k+1] * P
		F += Term
		normP = norm(P, Inf)
		P *= M
		small = eps()*norm(F, Inf)
		if norm(Term, Inf) <= small
			#=
			The termination condition is loosely inspired to the one in the paper.
			We check that the next lookahead terms in the series are small,
			estimating ||M^(k+r)|| with max(||M^k||, ||M^(k+1)||).
			Works well in practice, at least for exp, log, sqrt and pow.
			=#
			∆ = maximum(abs.(tay.coeffs[k+2:k+1+lookahead]))
			if ∆*max(normP, norm(P, Inf)) <= small
				break
			end
		end
	end
	return Matrix{C}(F)
end

function schurparlett(f::Func, A::Matrix{C}) where {Func, C<:Complex}
	if istriu(A)
		return schurparlett(f, A, eye(A))
	end

	T, Q, vals = schur(A)
	return schurparlett(f, T, Q)
end

function schurparlett(f::Func, T::Matrix{Comp}, Q::Matrix{Comp}) where {Func, Comp<:Complex}
	if isdiag(T)
		return Q*diagm(f.(diag(T)))*Q'
	end

	vals = diag(T)
	S, p = blockpattern(vals, 0.1)
	T, Q = copy(T), copy(Q)
	bsize = reorder!(T, Q, S, p)
	bend = cumsum(bsize)
	bbegin = bend - bsize .+ 1
	F = zeros(T)
	for j = 1:p
		J = bbegin[j]:bend[j]
		F[J,J] = atomicblock(f, T[J,J])
		for i = j-1:-1:1
			I = bbegin[i]:bend[i]
			K1, K2 = bbegin[i]:bend[j-1], bbegin[i+1]:bend[j]
			C = view(F,I,K1)*view(T,K1,J) - view(T,I,K2)*view(F,K2,J)
			Fij, scale = LAPACK.trsyl!('N', 'N', T[I,I], T[J,J], C, -1)
			F[I,J] = Fij/scale
		end
	end
	return Q*F*Q'
end
