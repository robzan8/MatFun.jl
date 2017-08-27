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
function reorder!(T::Matrix{C}, Q::Matrix{C}, S::Vector{Int64}, p::Int64) where {C <: Complex}
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
		minset = findmin(pos)[2]
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
