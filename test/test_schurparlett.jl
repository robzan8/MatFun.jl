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
	n = 1000
	srand(75)
	vals = zeros(n)
	vals[1] = n
	for i = 2:n
		vals[i] = vals[i-1] + rand()
	end
	vals *= 0.1
	T = Matrix{Complex128}(diagm(vals) + triu(randn(n, n), 1))
	@test cond(T) <= 1000

	@test MatFun.atomicblock(exp, T) ≈ LinAlg.expm(T)
	@test MatFun.atomicblock(log, T) ≈ LinAlg.logm(T)
	@test MatFun.atomicblock(sqrt, T) ≈ LinAlg.sqrtm(T)
	pow = (x) -> x^5.739
	@test MatFun.atomicblock(pow, T) ≈ pow(T)
end
