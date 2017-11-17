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
		V = svd(A, thin=true)[3]
		w = V[:, end]

		# Rational approximant on Z:
		num = C*(w.*f)
		den = C*w
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

	# Compute poles, residues and zeros:

	# Remove Froissart doublets:

	return
end

#=
Evaluate rational function in barycentric form.
=#
function reval(z::Vector{N}, f::Vector{N}, w::Vector{N}, x::N)::N where {
	N<:Union{Float32, Float64, Complex64, Complex128}}

	return reval(fill(x, 1))[1]
end
function reval(z::Vector{N}, f::Vector{N}, w::Vector{N}, x::A)::A where {
	N<:Union{Float32, Float64, Complex64, Complex128}, A<:Array{N}}
	v = x[:]
	C = 1./(v .- z.')     # Cauchy matrix
	r = (C*(w.*f))./(C*w) # result

	# Deal with input Inf as limit:
	r[isinf.(v)] = sum(w.*f)/sum(w)

	# Force interpolation at support points, to avoid Inf/Inf:
	for j = 1:length(v)
		i = findfirst(v[j] .== z)
		if i != 0
			r[j] = f[i]
		end
	end
	return reshape(r, size(x))
end
