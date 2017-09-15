#=
This file implements the algorithm described in:
"A Schur-Parlett Algorithm for Computing Matrix Functions"
(Davies and Higham, 2003)
=#

#=
blockpattern groups eigenvalues in sets, implementing Algorithm 4.1 of the paper.
vals is the vector of eigenvalues, delta is the tolerance.
The function returns S and p. S is the pattern, where S[i] = s means that
vals[i] has been assigned to set s. p is the number of sets identified.
=#
function blockpattern(vals::Vector{Num}, delta::Float64) where {Num<:Number}
	unassigned = -1
	S = fill(unassigned, length(vals))
	p = 0

	# assign vals[i] to set s.
	function assign(i::Int64, s::Int64)
		for k = i:length(vals)
			if vals[k] == vals[i]
				S[k] = s
			end
		end
	end
	# merge sets s and t.
	function mergesets(s::Int64, t::Int64)
		@assert s != t
		smin, smax = minmax(s, t)
		for k = 1:length(vals)
			if S[k] == smax
				S[k] = smin
			elseif S[k] > smax
				S[k] -= 1
			end
		end
		p -= 1
	end

	for i = 1:length(vals)
		if S[i] == unassigned
			p += 1
			assign(i, p)
		end
		for j = i+1:length(vals)
			if S[j] != S[i] && abs(vals[i]-vals[j]) <= delta
				if S[j] == unassigned
					assign(j, S[i])
				else
					mergesets(S[i], S[j])
				end
			end
		end
	end

	# complex conjugate eigenvalues can't be separated,
	# so we merge those sets:
	for i = 1:length(vals)-1
		j = i + 1
		if imag(vals[i]) != 0 && vals[i] == conj(vals[j]) && S[i] != S[j]
			mergesets(S[i], S[j])
		end
	end
	return S, p
end

#=
reorder reorders the complex Schur decomposition (T, Q) according to pattern S,
using the swapping strategy described in Algorithm 4.2.
The entries of S and vals are also reordered together with the corresponding eigenvalues.
The function returns vector blocksize: blocksize[i] is the size of the i-th
leading block on T's diagonal (after the reordering).
=#
function reorder!(T::Matrix{Num1}, Q::Matrix{Num1}, vals::Vector{Num2}, S::Vector{Int64}, p::Int64)
	where {Num1<:Number, Num2<:Number}
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

#=
parlettrec computes and returns f(T) using the Parlett recurrence.
T is block-triangular and blockend contains the indices at which each block ends.
The paper suggests to do the recurrence iteratively (Algorithm 5.1), we do it
recursively. This gives us better performance, as the recursion is cache oblivious.
Also, the iterative version suffered overheads in the non-uncommon case of many
small blocks.
=#
function parlettrec(f::Func, T::Matrix{Num}, blockend::Vector{Int64}) where {Func, Num<:Number}
	@assert length(blockend) > 0
	if length(blockend) == 1
		return atomicblock(f, T)
	end

	# split T in 2x2 superblocks of size ~n/2:
	n = size(T, 1)
	b = indmin(abs.(blockend .- n/2))
	bend = blockend[b]
	T11, T12 = T[1:bend, 1:bend], view(T, 1:bend, bend+1:n)
	T21, T22 = view(T, bend+1:n, 1:bend), T[bend+1:n, bend+1:n]

	F11 = parlettrec(f, T11, blockend[1:b])
	F22 = parlettrec(f, T22, blockend[b+1:end] .- bend)
	
	# Parlett says: T11*F12 - F12*T22 = F11*T12 - T12*F22
	F12, scale = LAPACK.trsyl!('N', 'N', T11, T22, F11*T12 - T12*F22, -1)

	return [F11 F12/scale; zeros(T21) F22]
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
