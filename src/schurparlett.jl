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
atomicblock computes f(T) using Taylor as described in Algorithm 2.6
=#
function atomicblock(f::Function, T::Matrix{Complex{R}}) where {R<:Real}
	n = size(T, 1)
	if n == 1
		return f.(T)
	end

	maxiter = 300
	μ = norm(UpperTriangular(eye(n) - abs.(triu(T, 1))) \ fill(1.0, n), Inf)
	σ = mean(diag(T))
	tay = f(σ + Taylor1(typeof(σ), max(10, n*3÷2)))
	M = UpperTriangular(T - diagm(fill(σ, n)))
	P = M
	F = UpperTriangular(diagm(fill(tay.coeffs[1], n)))
	for k = 1:maxiter
		needorder = k + n - 1
		@assert needorder <= tay.order + 1
		if needorder > tay.order
			tay = f(σ + Taylor1(typeof(σ), min(maxiter+n-1, tay.order*2)))
		end

		Term = tay.coeffs[k+1] * P
		F += Term
		small = eps()*norm(F, Inf)
		if norm(Term, Inf) <= small
			# we estimate ω[k+r] with |f⁽ᵏ⁺ʳ⁾(σ)| = |coeffs[k+p]|*Γ(k+p) where p = r + 1
			∆ = 0.0
			for p = 1:n
				∆ = max(∆, abs(tay.coeffs[k+p])*gamma(k+p)/gamma(p))
			end
			if μ*∆*norm(P, Inf) <= small
				break
			end
		end
		P *= M
	end
	return F
end
