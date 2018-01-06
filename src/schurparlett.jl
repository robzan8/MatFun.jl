#=
This file implements the algorithm described in:
"A Schur-Parlett Algorithm for Computing Matrix Functions"
(Davies and Higham, 2003)
(with some modifications to handle real matrices more efficiently).
=#

#=
blockpattern groups eigenvalues in sets, implementing Algorithm 4.1 of the paper.
The function returns S and p. S is the pattern, where S[i] = s means that
vals[i] has been assigned to set s. p is the number of sets identified.
=#
delta = 0.1
function blockpattern(vals::Vector{Complex128}, schurtype::Type)::Tuple{Vector{Int64}, Int64}
	unassigned = -1
	S = fill(unassigned, length(vals))
	p = 0

	# assign vals[i] to set s.
	function assign(i::Int64, s::Int64)
		lambda = vals[i]
		for k = i:length(vals)
			if vals[k] == lambda
				S[k] = s
			end
		end
	end
	function mergesets(s::Int64, t::Int64)
		if s == t
			return
		end
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
			if abs(vals[i]-vals[j]) <= delta
				if S[j] == unassigned
					assign(j, S[i])
				else
					mergesets(S[i], S[j])
				end
			end
		end
	end

	if schurtype <: Real
		# complex conjugate eigenvalues can't be separated in the
		# real Schur factorization, so we merge the sets involving them:
		i = 1
		while i < length(vals)
			if imag(vals[i]) != 0
				mergesets(S[i], S[i+1])
				i += 2
			else
				i += 1
			end
		end
	end
	return S, p
end

# analogous to ordschur/trsen:
function ordvec!(v::Vector{T}, select::B) where {T, B<:Union{Vector{Bool}, BitVector}}
	ilst = 1
	for ifst = 1:length(v)
		if select[ifst]
			if ifst != ilst
				v[ilst], v[ilst+1:ifst] = v[ifst], v[ilst:ifst-1]
			end
			ilst += 1
		end
	end
end

#=
reorder reorders the Schur decomposition (T, Q, vals) according to pattern S,
using the swapping strategy described in Algorithm 4.2.
The returned vector, blockend, contains the indices at which each block ends
(after the reordering).
=#
function reorder!(T::Matrix{N}, Q::Matrix{N}, vals::Vector{Complex{R}}, S::Vector{Int64}, p::Int64)::Vector{Int64} where {
	R<:Real, N<:Union{R, Complex{R}}}

	# for each set, calculate its mean position in S:
	pos = zeros(Float64, p)
	cou = zeros(Int64, p)
	for i = 1:length(S)
		set = S[i]
		pos[set] += i
		cou[set] += 1
	end
	pos ./= cou

	blockend = zeros(Int64, p)
	for set = 1:p
		numordered = (set == 1) ? 0 : blockend[set-1]
		minset = indmin(pos)
		select = [i <= numordered || S[i] == minset for i = 1:length(S)]
		ordschur!(T, Q, select)
		ordvec!(vals, select)
		ordvec!(S, select)
		blockend[set] = count(select)
		pos[minset] = Inf
	end
	return blockend
end

#=
taylorf computes f(T) using Taylor as described in Algorithm 2.6.
=#
function taylorf(f::Func, T::Mat, shift::N)::Mat where {
	Func, N<:Number, Mat<:Union{Matrix{N}, UpperTriangular{N, Matrix{N}}}}

	maxiter = 300
	lookahead = 10
	tay = f(shift + Taylor1(N, 20+lookahead))
	M = T - shift*Mat(eye(T))
	P = M
	F = tay.coeffs[1]*Mat(eye(T))
	for k = 1:maxiter
		needorder = k + lookahead
		@assert needorder <= tay.order + 1
		if needorder > tay.order
			tay = f(shift + Taylor1(N, min(maxiter+lookahead, tay.order*2)))
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
			Works well in practice, at least for exp, log and sqrt.
			=#
			delta = maximum(abs.(tay.coeffs[k+2:k+1+lookahead]))
			if delta*max(normP, norm(P, Inf)) <= small
				return F
			end
		end
	end
	error("Taylor did not converge.")
end

#=
blockf computes f(T) where T is an "atomic" block.
=#
function blockf(f::Func, T::Matrix{N}, vals::Vector{Complex{R}})::Matrix{N} where {
	Func, R<:Real, N<:Union{R, Complex{R}}}

	n = size(T, 1)
	if n == 1
		return f.(T)
	end
	if N <: Complex
		# T is complex and triangular with one cluster of eigenvalues.
		return taylorf(f, UpperTriangular(T), mean(vals))
	end
	if all(imag(vals) .== 0)
		# T is real and triangular with one cluster of real eigenvalues.
		return taylorf(f, UpperTriangular(T), mean(real(vals)))
	end
	if any(abs.(imag(vals)) .<= delta/2)
		# T is real and quasi-triangular with one cluster of complex
		# conjugate eigenvalues (and possibly some real eigenvalue).
		return taylorf(f, T, mean(real(vals)))
	end
	if n == 2
		#=
		T is real and 2x2 with two well-separated complex conjugate eigenvalues.
		We use the Hermite interpolating polynomial obtained with the
		Lagrange-Hermite formula, simplified assuming f(conj(x)) == conj(f(x)).
		=#
		v1, v2 = vals[1], vals[2]
		@assert v1 == conj(v2)
		psi1 = f(v1)/(v1 - v2)
		return 2.0*real(psi1*(T - [v2 0; 0 v2]))
	end
	#=
	T is real with two well-separated, complex conjugate clusters of eigenvalues.
	We use the complex Schur factorization to further divide T.
	Doing complex Schur at the block level is more efficient than using complex
	arithmetic for the whole matrix.
	=#
	U, Z, vals = schur(Matrix{Complex{R}}(T))::Tuple{Matrix{Complex{R}}, Matrix{Complex{R}}, Vector{Complex{R}}}
	select = imag(vals) .> 0
	ordschur!(U, Z, select)
	ordvec!(vals, select)
	@assert count(select) == n÷2
	F = Z*recf(f, U, vals, [n÷2, n])*Z'
	realF = real(F)
	if norm(imag(F), Inf) > sqrt(eps())*norm(realF, Inf)
		warn("T is real but f(T) has non-negligible imaginary part.")
		# note: we need f(T) to be real, because trsyl can't handle
		# quasi-triangular complex matrices.
	end
	return realF
end

#=
recf computes f(T) using the Parlett recurrence (T is block-triangular).
The paper suggests to do the recurrence iteratively (Algorithm 5.1),
we do it recursively as this version performs better.
=#
function recf(f::Func, T::Matrix{N}, vals::Vector{Complex{R}}, blockend::Vector{Int64})::Matrix{N} where {
	Func, R<:Real, N<:Union{R, Complex{R}}}

	@assert length(blockend) > 0
	if length(blockend) == 1
		return blockf(f, T, vals)
	end

	# split T in 2x2 superblocks of size ~n/2:
	n = size(T, 1)
	b = indmin(abs.(blockend .- n/2))
	bend = blockend[b]
	T11, T12 = T[1:bend, 1:bend], view(T, 1:bend, bend+1:n)
	T21, T22 = view(T, bend+1:n, 1:bend), T[bend+1:n, bend+1:n]

	F11 = recf(f, T11, vals[1:bend], blockend[1:b])
	F22 = recf(f, T22, vals[bend+1:n], blockend[b+1:end] .- bend)
	
	# T11*F12 - F12*T22 = F11*T12 - T12*F22
	F12, scale = LAPACK.trsyl!('N', 'N', T11, T22, F11*T12 - T12*F22, -1)

	return [F11 F12/scale; zeros(T21) F22]
end

#=
Computes f(A) using the Schur-Parlett algorithm.
When A is a real matrix, computation will be done mostly in real arithmetic
and the algorithm will assume f(conj(x)) == conj(f(x)).
=#
function schurparlett(f::Func, A::Matrix{N})::Matrix{N} where {Func, N<:Union{Float64, Complex128}}
	T, Q, vals = schur(A)
	return schurparlett(f, T, Q, vals)
end

#=
Analogous to schurparlett(f, A), but A is provided as its Schur decomposition.
=#
function schurparlett(f::Func, T::Matrix{Float64}, Q::Matrix{Float64}, vals::Vector{Float64})::Matrix{Float64} where {Func}
	return schurparlett(f, T, Q, Vector{Complex128}(vals))
end
function schurparlett(f::Func, T::Matrix{N}, Q::Matrix{N}, vals::Vector{Complex128})::Matrix{N} where {
	Func, N<:Union{Float64, Complex128}}

	if size(T, 1) != size(T, 2)
		throw(DimensionMismatch("T must be square"))
	end
	if size(Q) != size(T) || length(vals) != size(T, 1)
		throw(DimensionMismatch("T, Q, vals dimension mismatch"))
	end
	if size(T, 1) == 0
		error("T is empty")
	end

	d = diag(T)
	D = diagm(d)
	if norm(T-D, Inf) <= 1000*eps()*norm(D, Inf) # T is diagonal
		return Q*diagm(f.(d))*Q'
	end

	S, p = blockpattern(vals, N)
	blockend = reorder!(T, Q, vals, S, p)
	return Q*recf(f, T, vals, blockend)*Q'
end
