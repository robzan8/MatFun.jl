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
function poles_to_moebius!(p::Vector{Complex{R}})::Tuple{Vector{R}, Vector{R}, Vector{R}} where {R<:Real}
	mu = ones(length(p))
	rho = ones(length(p))
	eta = zeros(length(p))

	sel = p .== 0

	p[sel] = 1
	mu[sel] = Inf
	rho[sel] = 0
	eta[sel] = 1

	return mu, rho, eta
end

function ratkrylov(A::Matrix{N}, b::Vector{N}, p::Vector{C}) where {N<:Number, C<:Complex}
	B = eye(A)
	m, n = length(p), length(b)
	V = zeros(N, n, m+1)
	K = zeros(N, m+1, m)
	H = zeros(K)
	realopt = N <: Real

	# Cannot use real arithmetic if the poles are not ordered canonically.
	if realopt && !canonical_cplx(p)
		realopt = false
		println("can't use real arithmetic")
	end

	mu, rho, eta = poles_to_moebius!(p)

	# run_krylov:
	bd = false
	bd_tol = eps()
	# Starting vector.
	V[:, 1] = b/norm(b)
	j = 1
	while j <= m
		# Computing the continuation combination.
		if j == 1
			U[1, 1] = 1
		else
			Q = qr(K[1:j, 1:j-1]/mu[j] - H[1:j, 1:j-1]/p[j], thin=false)[1]
			U[1:j, j] = Q[:, end]
		end

		# Compute new vector.
		w = V[:, 1:j]*U[1:j, j]
		w = rho[j]*(A*w) - eta[j]*w
		w = (B/mu[j] - A/p[j]) \ w

		# Orthogonalization.
		if realopt && !isreal(p[j])
			println("yay")
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
			u, h = U[1:j-1, j-1], H[1:j+1, j-1:j]
			rcnt = [real(u) imag(u); 0 0; 0 0]
			K[1:j+1, j-1:j] = h*cp + rho[j-1]*rcnt
			H[1:j+1, j-1:j] = h/mu[j-1] + eta[j-1]*rcnt
		else
			#w = (B/mu[j] - A/N(p[j])) \ w
			# MGS
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
			u, h = [U[1:j, j]; 0], H[1:j+1, j]
			println(u)
			K[1:j+1, j] = rho[j]*u + h/p[j]
			H[1:j+1, j] = eta[j]*u + h/mu[j]
		end # realopt

		j = j+1
	end # while j <= m

	if bd == true
		println("lucky breakdown at iteration ", j)
		V = V[:, 1:j]
		K = K[1:j, 1:j-1]
		H = H[1:j, 1:j-1]
	end
	return V, K, H
end
