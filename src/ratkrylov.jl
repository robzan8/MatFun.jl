#=
canonical_cplx checks if the poles p are ordered canonically.
=#
function canonical_cplx(p::Vector{C})::Bool where {C<:Complex}
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
function poles_to_moebius(p::Vector{Complex{R}})::
	Tuple{Vector{Complex{R}}, Vector{R}, Vector{R}, Vector{R}} where {R<:Real}
	
	p = copy(p)
	mu = ones(R, length(p))
	rho = ones(R, length(p))
	eta = zeros(R, length(p))

	sel = p .== 0

	p[sel] = 1
	mu[sel] = Inf
	rho[sel] = 0
	eta[sel] = 1

	return p, mu, rho, eta
end

function ratkrylov(A::Mat, b::Vector{N}, p::Vector{Complex{R}}) where {
	R<:Union{Float32, Float64}, N<:Union{R, Complex{R}}, Mat<:Union{Matrix{N}, SparseMatrixCSC{N}}}

	B = eye(A)
	m, n = length(p), length(b)
	V = zeros(N, n, m+1)
	K = zeros(N, m+1, m)
	H = zeros(K)
	realopt = N <: Real

	# Cannot use real arithmetic if the poles are not ordered canonically.
	if realopt && !canonical_cplx(p)
		error("can't use real arithmetic")
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
			u = ones(Complex{R}, 1) # U[1:j, j] in the original code
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
						v = V[:, reo_i]
						hh = V[:, j+1]'*v
						V[:, j+1] -= v*hh
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
					v = V[:, reo_i]
					hh = w'*v
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
		V = V[:, 1:j]
		K = K[1:j, 1:j-1]
		H = H[1:j, 1:j-1]
	end
	return V, K, H
end
