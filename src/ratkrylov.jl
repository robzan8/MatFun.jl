#=
canonical_cplx checks if the poles p are ordered canonically.
=#
function canonical_cplx(p::Vector{Complex128})
	m = length(p)
	j = 1
	while j <= m
		if isreal(p[j]) || isinf(p[j])
			j = j+1
		else
			if j == m || p[j+1] != conj(p[j])
				return false
			end
			j = j+2
		end
	end
	return true
end

#=
Moebius transformation with poles p.
p != 0 is replaced by (p, mu) := (p,   1) and (rho, eta) := (1, 0),
p  = 0 is replaced by (p, mu) := (1, Inf) and (rho, eta) := (0, 1).
=#
function poles_to_moebius(p::Vector{Complex128})
	p = copy(p)
	mu = ones(length(p))
	rho = ones(length(p))
	eta = zeros(length(p))

	sel = p .== 0

	p[sel] = 1
	mu[sel] = Inf
	rho[sel] = 0
	eta[sel] = 1

	return p, mu, rho, eta
end

function ratkrylov(A::Mat, b::Vector{N}, p::Vector{Complex128}) where {
	N<:Union{Float64, Complex128}, Mat<:Union{Matrix{N}, SparseMatrixCSC{N}}}

	if size(A, 1) != size(A, 2)
		throw(DimensionMismatch("A must be square"))
	end
	if size(A, 2) != length(b)
		throw(DimensionMismatch("can't do A*b"))
	end
	m, n = length(p), length(b)
	if n == 0
		error("A is empty")
	end
	if m == 0
		error("no poles")
	end
	if m > n
		error("too many poles")
	end

	B = eye(A)
	V = zeros(N, n, m+1)
	K = zeros(N, m+1, m)
	H = zeros(K)
	realopt = N <: Real

	if realopt && !canonical_cplx(p)
		error("Cannot use real arithmetic if the poles are not ordered canonically.")
	end

	p, mu, rho, eta = poles_to_moebius(p)

	# run_krylov:
	bd = false
	bd_tol = eps()
	# Starting vector.
	V[:, 1] = b/norm(b)
	j = 1
	while j <= m
		if realopt && !isreal(p[j])
			# Computing the continuation combination.
			u = ones(Complex128, 1) # U[1:j, j] in the original code
			if j > 1
				Q = qr(K[1:j, 1:j-1]/mu[j] - H[1:j, 1:j-1]/p[j], thin=false)[1]
				u = Q[:, end]
			end

			# Compute new vector.
			w = V[:, 1:j]*u
			w = rho[j]*(A*w) - eta[j]*w
			w = (B/mu[j] - A/p[j]) \ w

			# Orthogonalization.
			V[:, j+1] = real(w)
			V[:, j+2] = imag(w)
			# MGS
			for j = j:j+1
				for reo = 0:1
					for reo_i = 1:j
						vi = view(V, 1:size(V, 1), reo_i)
						vj = view(V, 1:size(V, 1), j+1)
						hh = vi'*vj
						V[:, j+1] = vj - vi*hh
						H[reo_i, j] += hh
					end
				end
				normw = norm(V[:, j+1])
				H[j+1, j] = normw
				V[:, j+1] /= normw
				if normw < bd_tol*norm(H[1:j, j])
					bd = true
					break
				end
			end

			# Setting the decomposition.
			rp, ip = real(1/p[j-1]), imag(1/p[j-1])
			cp = [rp ip; -ip rp]
			u0, h = [real(u) imag(u); 0 0; 0 0], H[1:j+1, j-1:j]
			K[1:j+1, j-1:j] = rho[j-1]*u0 + h*cp
			H[1:j+1, j-1:j] = eta[j-1]*u0 + h/mu[j-1]
		else
			pj = N(p[j])
			# Computing the continuation combination.
			u = ones(N, 1) # U[1:j, j] in the original code
			if j > 1
				Q = qr(K[1:j, 1:j-1]/mu[j] - H[1:j, 1:j-1]/pj, thin=false)[1]
				u = Q[:, end]
			end

			# Compute new vector.
			w = V[:, 1:j]*u
			w = rho[j]*(A*w) - eta[j]*w
			w = (B/mu[j] - A/pj) \ w

			# Orthogonalization, MGS.
			for reo = 0:1
				for reo_i = 1:j
					v = view(V, 1:size(V, 1), reo_i)
					hh = v'*w
					w -= v*hh
					H[reo_i, j] += hh
				end
			end
			normw = norm(w)
			H[j+1, j] = normw
			V[:, j+1] = w/normw
			if normw < bd_tol*norm(H[1:j, j])
				bd = true
				break
			end

			# Setting the decomposition.
			u0, h = [u; 0], H[1:j+1, j]
			K[1:j+1, j] = rho[j]*u0 + h/pj
			H[1:j+1, j] = eta[j]*u0 + h/mu[j]
		end # realopt

		j = j+1
	end # while j <= m

	if bd == true
		println("lucky breakdown")
		V = V[:, 1:j]
		K = K[1:j, 1:j-1]
		H = H[1:j, 1:j-1]
	end
	return V, K, H
end

function ratkrylovf(f::Func, A::Mat, b::Vector{N}, p::Vector{Complex128}) where {
	Func, N<:Union{Float64, Complex128}, Mat<:Union{Matrix{N}, SparseMatrixCSC{N}}}

	V, K, H = ratkrylov(A, b, p)

	m = size(V, 2) - 1 # may be < length(p) in case of breakdown
	Am = isinf(p[m]) ? [H/K[1:m,1:m] V'*(A*V[:,end])] : V'*A*V

	return V*(schurparlett(f, Am)*(V'*b))
end

#=
automatic poles selection:
=#
function ratkrylovf(f::Func, A::Mat, b::Vector{N}, mmax::Int64=100, tol::Float64=1e-13) where {
	Func, N<:Union{Float64, Complex128}, Mat<:Union{Matrix{N}, SparseMatrixCSC{N}}}

	rad = Mat<:SparseMatrixCSC ? min(norm(A, 1), norm(A, Inf), vecnorm(A)) : norm(A, 2)
	rad = max(rad, sqrt(eps()))

	nsamples = Mat<:SparseMatrixCSC ? nnz(A)÷2 : prod(size(A))
	nsamples = max(nsamples, 100)

	pol = aaa(f, lpdisk(rad, nsamples), tol, mmax)[2]
	if N <: Real
		# Assume f(conj(x)) == conj(f(x)),
		# real ratkrylov wants conjugated and canonically ordered poles.
		pos = pol[imag.(pol) .> 1e-13]
		pol = pol[abs.(imag.(pol)) .<= 1e-13]
		pol[1:end] = real.(pol)
		for i = 1:length(pos)
			pol = [pol; pos[i]; conj(pos[i])]
		end
	end
	# Inf as last pole allows for faster projection of A in the Krylov space:
	pol = [pol; Inf]

	return ratkrylovf(f, A, b, pol)
end
