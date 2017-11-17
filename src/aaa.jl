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
		A = F.*C - C.*transpose(f)
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
end
