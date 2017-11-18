using MatFun, Base.Test

@testset "aaa" begin
	
end

@testset "poles_to_moebius" begin
	p, mu, rho, eta = MatFun.poles_to_moebius([1.0+0im, 0, Inf*im, 0, 47+79im])
	@test p == [1, 1, Inf*im, 1, 47+79im]
	@test mu == [1, Inf, 1, Inf, 1]
	@test rho == [1, 0, 1, 0, 1]
	@test eta == [0, 1, 0, 1, 0]
end

@testset "ratkrylov" begin
	srand(789)
	n = 10
	m = 5
	for sparse = [false, true]
		A = randn(n, n)
		if sparse
			A = SparseMatrixCSC(A)
		end
		b = randn(n)
		p = randn(m) + im*randn(m)
		@test_throws ErrorException ratkrylov(A, b, p) # poles not canonically ordered
		p = [1.0+0im, Inf, 3+4im, 3-4im, 5]
		V, K, H = ratkrylov(A, b, p)
		@test A*V*K ≈ V*H
		p = [1.0+0im, 2, 3+4im, 3-4im, Inf]
		V, K, H = ratkrylov(A, b, p)
		@test A*V[:,1:m]*K[1:m,:] ≈ V*H
		Am = H/K
		Am[:,end] = V'*(A*V[:,end])
		@test V'*A*V ≈ Am
		p = fill(Inf+0im, m)
		V, K, H = ratkrylov(A, b, p)
		@test A*V[:,1:m]*K[1:m,:] ≈ V*H && K[1:m,:] == eye(m)

		A = randn(n, n) + im*randn(n, n)
		if sparse
			A = SparseMatrixCSC(A)
		end
		b = randn(n) + im*randn(n)
		p = randn(m) + im*randn(m)
		V, K, H = ratkrylov(A, b, p)
		@test A*V*K ≈ V*H
		p[m] = Inf
		V, K, H = ratkrylov(A, b, p)
		@test A*V[:,1:m]*K[1:m,:] ≈ V*H
		Am = H/K
		Am[:,end] = V'*(A*V[:,end])
		@test V'*A*V ≈ Am #it's failing!!!!!!!!!!!!!!!!!!!!
		p = fill(Inf+0im, m)
		V, K, H = ratkrylov(A, b, p)
		@test A*V[:,1:m]*K[1:m,:] ≈ V*H && K[1:m,:] == eye(m)
	end
end
