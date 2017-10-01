using MatFun, Base.Test

@testset "blockpattern" begin
	vals = zeros(Complex128, 10)
	S, p = MatFun.blockpattern(vals, Complex128)
	@test S == fill(1, length(vals))
	@test p == 1

	vals = Vector{Complex128}(collect(1:10))
	S, p = MatFun.blockpattern(vals, Complex128)
	@test S == vals
	@test p == length(vals)

	vals = Vector{Complex128}(randperm(20)*0.09)
	S, p = MatFun.blockpattern(vals, Complex128)
	@test S == fill(1, length(vals))
	@test p == 1

	vals = [1, 1+0.09im, 2, 3, 1+0.19im, 2+0.09im, 2+0.19im, 3, 3, 3+0.09im, 3+0.19im, 1, 1+0.29im]
	S, p = MatFun.blockpattern(vals, Complex128)
	@test S == real.(vals)
	@test p == 3

	vals = [1, 1+4im, 1-4im, 3, 1+4.05im, 1-4.05im, 2, 3, 3, 3+0.09im, 3-0.09im, 1, 1]
	S, p = MatFun.blockpattern(vals, Float64)
	@test S[1] != S[2] && S[2] == S[3] && S[2] == S[5] && S[5] == S[6]
	@test S[10] == S[11] && S[10] == S[9] && S[10] != S[12]
	@test p == 4
end

@testset "reorder" begin
	vals = Vector{Complex128}([1, 1, 1, 1, 1, 2, 1, 2, 3, 2, 2, 3, 3, 1, 3])
	T = diagm(vals)
	Q = eye(T)
	A = copy(T)
	S, p = MatFun.blockpattern(vals, Complex128)
	blockend = MatFun.reorder!(T, Q, vals, S, p)
	newvals = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
	@test T == diagm(newvals)
	@test vals == newvals
	@test Q*T*Q' == A
	@test S == newvals
	@test blockend == cumsum([7, 4, 4])

	srand(42)
	A = Matrix{Complex128}(randn(20, 20)/18)
	@test cond(A) <= 1000
	T, Q, vals = schur(A)
	S, p = MatFun.blockpattern(vals, Complex128)
	MatFun.reorder!(T, Q, vals, S, p)
	@test vals == diag(T)
	for i = 1:length(vals)÷2
		j = length(vals) - i + 1
		if S[i] != S[j]
			@test abs(vals[i]-vals[j]) > MatFun.delta
		end
	end
	@test Q*T*Q' ≈ A
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
