#=
comment
=#
function aaa(func::Func, Z::Vector{N}, tol::Float64=1e-13, mmax::Int64=100) where {
	Func, N<:Union{Float32, Float64, Complex64, Complex128}}

	F = (Func <: Vector) ? Vector{N}(func) : Vector{N}(func.(Z))

	# Remove any infinite or NaN function values (avoid SVD failures):
	keep = isfinite.(F)
	Z = Z[keep]
	F = F[keep]
	M = length(Z)
	@assert M == length(F)

	abstol = tol*norm(F, Inf)
	J = collect(1:M)            # indices of the non-support points
	z = Vector{N}(0)            # support points
	f = Vector{N}(0)            # corresponding data values
	w = Vector{N}(0)            # weights
	C = Matrix{N}(M, 0)         # Cauchy matrix
	errvec = real(Vector{N}(0))
	R = fill(mean(F), M)

	for m = 1:mmax
		# Select next support point where error is largest:
		j = indmax(abs.(F - R))
		z = [z; Z[j]]
		f = [f; F[j]]
		J = J[J .!= j]
		C = [C 1./(Z .- Z[j])]

		# Compute weights:
		A = F.*C - C.*f.'
		V = svd(A[J,:], thin=true)[3]
		w = V[:,m]

		# Rational approximant on Z:
		num, den = C*(w.*f), C*w
		R = copy(F)
		R[J] = num[J]./den[J]

		err = norm(F - R, Inf)
		errvec = [errvec; err]
		if err <= abstol
			break
		end
	end

	# Remove support points with zero weight:
	keep = w .!= 0
	z = z[keep]
	f = f[keep]
	w = w[keep]

	# Construct function handle:
	r = (x) -> reval(z, f, w, x)

	# Compute poles and zeros via generalized eigenvalues:
	m = length(w)
	B = eye(N, m+1)
	B[1,1] = 0
	pol = eigvals([0 w.'; ones(N, m, 1) diagm(z)], B)
	zer = eigvals([0 (w.*f).'; ones(N, m, 1) diagm(z)], B)
	#=
	Note: some Inf poles can come up as NaN,
	since in Julia (1.+2im)/(0.+0im) is NaN+NaN*im
	(https://github.com/JuliaLang/julia/issues/9790).
	=#
	pol = pol[isfinite.(pol)]
	zer = zer[isfinite.(zer)]
	@assert length(pol) == m-1 && length(zer) == m-1

	# Compute residues via discretized Cauchy integral:
	dz = (1e-5)*exp.(2im*pi*collect(1:4)/4)
	res = r(pol .+ dz.')*(dz/4)

	# We don't remove numerical Froissart doublets,
	# which are rare anyway if aaa is used correctly. Sorry.

	return r, pol, res, zer, z, f, w, errvec
end

#=
Evaluate rational function in barycentric form.
=#
function reval(z::Vector{N}, f::Vector{N}, w::Vector{N}, x::Vector{X}) where {N<:Number, X<:Number}
	C = 1./(x .- z.')     # Cauchy matrix
	r = (C*(w.*f))./(C*w) # result

	# Deal with input Inf as limit:
	r[isinf.(x)] = sum(w.*f)/sum(w)

	# Force interpolation at support points, to avoid Inf/Inf:
	for j = 1:length(x)
		i = findfirst(x[j] .== z)
		if i != 0
			r[j] = f[i]
		end
	end
	return r
end
function reval(z::Vector{N}, f::Vector{N}, w::Vector{N}, x::X) where {N<:Number, X<:Number}
	return reval(z, f, w, [x])[1]
end
function reval(z::Vector{N}, f::Vector{N}, w::Vector{N}, A::Array{X}) where {N<:Number, X<:Number}
	return reshape(reval(z, f, w, A[:]), size(A))
end

#=
Evaluate rational matrix function in barycentric form.
For experiments only, use schurparlett instead.
Does not work when A has some support point as eigenvalue.
=#
function revalm(z::Vector{N}, f::Vector{N}, w::Vector{N}, A::Matrix{X}) where {N<:Number, X<:Number}
	Num = zeros(A)
	Den = zeros(A)
	for j = 1:length(z)
		W = (w[j]*eye(A)) / (A - z[j]*eye(A))
		Num += W*f[j]
		Den += W
	end
	return Num/Den
end

function lppoints(k::Int64)::Vector{Complex128}
	N = 1 << k
	p = collect(0:N-1)/(N+0im)
	for i = 0:N-1
		for j = 0:k-1
			p[i+1] += im*(count_ones(i >> j) & 1)/(2 << j)
		end
	end
	return p
end

function lpdisk(rad::R, n::Int64)::Vector{Complex{R}} where {R<:Real}
	@assert n >= 8
	k = Int64(ceil(log2(n*4/pi)))
	p = lppoints(k)
	N = 1 << k
	scale = sqrt(N*rad*rad*pi/n)
	p = (p .- (0.5+0.5im))*scale
	return p[abs2.(p) .<= rad*rad]
end
