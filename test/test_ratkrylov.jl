using MatFun, Base.Test

@testset "canonical_cplx" begin
	@test MatFun.canonical_cplx([1.0+0im, 2.0, 3.0]) == true
	@test MatFun.canonical_cplx([1.0+7im, 1.0-7im, 3.0]) == true
	@test MatFun.canonical_cplx([1.0+6im, 1.0+3im, 3.0]) == false
	@test MatFun.canonical_cplx([1.0+0im, 2.0, 3.0+3im]) == false
end

@testset "poles_to_moebius" begin
	p
end

@testset "taylorf" begin
	for t = 1:2
		n = 10
		tol = 100*eps()
		srand(75*t)
		vals = zeros(Complex128, n)
		vals[1] = 50
		for i = 2:n
			vals[i] = (t == 1) ? vals[1]+exp(im*2*pi*i/n) : vals[i-1]+rand()
		end
		vals *= 0.1
		shift = mean(vals)
		T = diagm(vals) + triu(randn(n, n), 1)
		@test cond(T) <= 1000

		F1 = LinAlg.expm(T)
		F2 = (t == 1)?
			MatFun.taylorf(exp, UpperTriangular(T), shift):
			MatFun.taylorf(exp, Matrix{Float64}(T), Float64(shift))
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.logm(T)
		F2 = (t == 1)?
			MatFun.taylorf(log, UpperTriangular(T), shift):
			MatFun.taylorf(log, Matrix{Float64}(T), Float64(shift))
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.sqrtm(T)
		F2 = (t == 1)?
			MatFun.taylorf(sqrt, UpperTriangular(T), shift):
			MatFun.taylorf(sqrt, Matrix{Float64}(T), Float64(shift))
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol
	end
end

@testset "blockf" begin
	srand(911)

	vals = [7.8 + 0im + rand()]
	T = eye(1)*vals[1]
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, Matrix{Float64}(T), vals)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals)

	vals = [5.0+0im, 5.1, 5.2, 4.9]
	T = diagm(vals) + im*triu(randn(4, 4), 1)
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, T, vals)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals)

	vals = [5.0, 5.1, 5.2, 4.9]
	T = diagm(vals) + triu(randn(4, 4), 1)
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, T, vals.+0im)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals.+0im)

	vals = [5.0, 5.1, 5.2, 4.9] .+ im*MatFun.delta*0.4
	T = diagm(real(vals)) + triu(randn(4, 4), 1)
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, T, vals)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals)

	T = [7.0 4.0+rand(); -2.0+rand() 7.0]
	U, Z, vals = schur(T)
	@test vals[1] == conj(vals[2]) && all(abs.(imag(vals)) .> MatFun.delta*0.5)
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, T, vals)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals)

	T2 = T + randn(2, 2)*MatFun.delta*0.1
	T = [T randn(2, 2); zeros(T) T2]
	U, Z, vals = schur(T)
	@test vals[1] == conj(vals[2]) && vals[3] == conj(vals[4]) && all(abs.(imag(vals)) .> MatFun.delta*0.5)
	@test LinAlg.expm(T) ≈ MatFun.blockf(exp, T, vals)
	@test LinAlg.sqrtm(T) ≈ MatFun.blockf(sqrt, T, vals)
end

@testset "schurparlett" begin
	for t = 1:3
		n = 10
		tol = 100*eps()
		srand(666*t)
		A::Union{Matrix{Float64}, Matrix{Complex128}} = Matrix{Float64}(0, 0)
		if t == 1
			A = diagm(randn(n) .+ 5)
		elseif t == 2
			A = 5*eye(n) + randn(n, n)/30 + im*randn(n, n)/30
		else
			A = 5*eye(n) + randn(n, n)/15
		end
		@test cond(A) <= 1000

		F1 = LinAlg.expm(A)
		F2 = schurparlett(exp, A)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.logm(A)
		F2 = schurparlett(log, A)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.sqrtm(A)
		F2 = schurparlett(sqrt, A)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol
	end
end
