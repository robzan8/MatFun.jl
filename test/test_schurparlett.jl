using MatFun, Base.Test

@testset "blockpattern" begin
	vals = zeros(Complex128, 10)
	S, p = MatFun.blockpattern(vals, 0.1)
	@test S == fill(1, length(vals))
	@test p == 1

	vals = Vector{Complex128}(collect(1:10))
	S, p = MatFun.blockpattern(vals, 0.1)
	@test S == vals
	@test p == length(vals)

	vals = [1, 1+0.1im, 2, 3, 1+0.2im, 2+0.1im, 2+0.2im, 3, 3, 3+0.1im, 3+0.2im, 1, 1+0.3im]
	S, p = MatFun.blockpattern(vals, 0.10001)
	@test S == real.(vals)
	@test p == 3

	vals = Vector{Complex128}(randperm(20))
	S, p = MatFun.blockpattern(vals, 1.0)
	@test S == fill(1, length(vals))
	@test p == 1
end

@testset "reorder" begin
	vals = Vector{Complex128}([1, 1, 1, 1, 1, 2, 1, 2, 3, 2, 2, 3, 3, 1, 3])
	T = diagm(vals)
	Q = eye(Complex128, length(vals))
	A = copy(T)
	S, p = MatFun.blockpattern(vals, 0.1)
	blocksize = MatFun.reorder!(T, Q, S, p)
	newvals = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
	@test T == diagm(newvals)
	@test Q*T*Q' == A
	@test S == newvals
	@test blocksize == [7, 4, 4]

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
	@test Q*T*Q' ≈ A
end

@testset "atomicblock" begin
	for t = 1:2
		n = 10
		tol = 30
		srand(75*t)
		vals = zeros(Complex128, n)
		vals[1] = n
		for i = 2:n
			vals[i] = (t == 1) ? vals[1]+exp(im*2*pi*i/n) : vals[i-1]+rand()
		end
		vals *= 0.1
		T = diagm(vals) + triu(randn(n, n), 1)
		@test cond(T) <= 1000

		F1 = LinAlg.expm(T)
		F2 = MatFun.atomicblock(exp, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol*eps()

		F1 = LinAlg.logm(T)
		F2 = MatFun.atomicblock(log, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol*eps()

		F1 = LinAlg.sqrtm(T)
		F2 = MatFun.atomicblock(sqrt, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol*eps()

		pow = (x) -> x^Float64(pi)
		F1 = pow(T)
		F2 = MatFun.atomicblock(pow, T)
		relerr = norm(F2-F1)/norm(F1)
		@test relerr <= tol*eps()
	end
end
