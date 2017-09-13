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

function parlettrec(f::Func, T::Matrix{Num}, blockend::Vector{Int64}) where {Func, Num<:Number}
	@assert length(blockend) > 0
	if length(blockend) == 1
		return atomicblock(f, T)
	end

	# Split T in 2x2 superblocks of size ~n/2:
	n = size(T, 1)
	b = indmin(abs.(blockend .- n/2))
	bend = blockend[b]
	T11, T12 = T[1:bend,     1:bend], T[1:bend,     bend+1:end]
	T21, T22 = T[bend+1:end, 1:bend], T[bend+1:end, bend+1:end]

	F11 = parlettrec(f, T11, blockend[1:b])
	F22 = parlettrec(f, T22, blockend[b+1:end] .- bend)
	
	# Parlett says: T11*F12 - F12*T22 = F11*T12 - T12*F22
	F12, scale = LAPACK.trsyl!('N', 'N', T11, T22, F11*T12 - T12*F22, -1)

	return [F11 F12/scale; T21#=0=# F22]
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
	blocksize = reorder!(T, Q, S, p)
	F = parlettrec(f, T, cumsum(blocksize))
	return Q*F*Q'
end
