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
	#@test Q*T*Q' == A
	#@test S == newvals
	#@test blocksize == cumsum([7, 4, 4])
#=
	srand(42)
	A = Matrix{Complex128}(randn(20, 20))
	@test cond(A) <= 1000
	T, Q, vals = schur(A)
	delta = 1.8
	S, p = MatFun.blockpattern(vals, delta)
	blocksize = MatFun.reorder!(T, Q, S, p)
	newvals = diag(T)
	a = 1
	for set = 1:p
		b = a + blocksize[set] - 1
		for i = a+1:b
			@test S[a] == S[i] && abs(newvals[a]-newvals[i]) <= delta*blocksize[set]+sqrt(eps())
		end
		a = b + 1
	end
	for i = 1:length(newvals)÷2
		j = length(newvals) - i + 1
		if S[i] != S[j]
			@test abs(newvals[i]-newvals[j]) > delta
		end
	end
	@test Q*T*Q' ≈ A=#
end
#=
@testset "atomicblock" begin
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
		T = diagm(vals) + triu(randn(n, n), 1)
		@test cond(T) <= 1000

		F1 = LinAlg.expm(T)
		F2 = MatFun.atomicblock(exp, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.logm(T)
		F2 = MatFun.atomicblock(log, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		F1 = LinAlg.sqrtm(T)
		F2 = MatFun.atomicblock(sqrt, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol

		pow = (x) -> x^Float64(pi)
		F1 = pow(T)
		F2 = MatFun.atomicblock(pow, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol
	end
end

@testset "schurparlett" begin
	for t = 1:3
		n = 10
		tol = 100*eps()
		srand(666*t)
		A = Matrix{Complex128}(0, 0)
		if t == 1
			A = 5*eye(n) + randn(n, n) + im*randn(n, n)
		elseif t == 2
			A = Matrix{Complex128}(diagm(randn(n) .+ 5))
		else
			vals = [4.0, 4, 4, 5, 5, 5, 7, 7, 7, 7]
			for i = 1:length(vals)
				vals[i] += rand()*0.1
			end
			A = Matrix{Complex128}(diagm(vals) + triu(randn(n, n), 1))
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

		pow = (x) -> x^Float64(pi)
		F1 = pow(A)
		F2 = schurparlett(pow, A)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol
	end
end
=#