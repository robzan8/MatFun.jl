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
	S, p = MatFun.blockpattern(vals, 0.1)
	@test S == real.(vals)
	@test p == 3

	vals = Vector{Complex128}(randperm(20))
	S, p = MatFun.blockpattern(vals, 1.0)
	@test S == fill(1, length(vals))
	@test p == 1
end
