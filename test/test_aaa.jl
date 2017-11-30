using MatFun, Base.Test

@testset "aaa" begin
	Z = collect(linspace(-3.1, 3.1, 200))
	r, pol, res, zer, z, f, w, errvec = aaa(gamma, Z)
	@test r(Z) ≈ gamma.(Z)
	@test any(pol .≈ -1) && any(pol .≈ -2) && any(pol .≈ -3)
	@test norm(gamma.(zer), Inf) <= 0.1

	srand(789)
	poles = randn(10) + im*randn(10)
	residues = randn(10) + im*randn(10)
	func = (x) -> sum(residues./(x .- poles))
	Z = MatFun.lpdisk(sqrt(2.0), 400)
	F = func.(Z)
	r, pol, res, zer, z, f, w, errvec = aaa(F, Z)
	@test r(Z) ≈ F
	@test sort(real(pol)) ≈ sort(real(poles)) && sort(imag(pol)) ≈ sort(imag(poles))
	@test sort(real(res)) ≈ sort(real(residues)) && sort(imag(res)) ≈ sort(imag(residues))
	@test norm(func.(zer), Inf) <= 0.1
end

# testset matrix reval
