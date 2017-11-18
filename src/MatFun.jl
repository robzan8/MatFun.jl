__precompile__()

module MatFun

using Base.LinAlg, TaylorSeries

export schurparlett, aaa, ratkrylov, ratkrylovf

include("schurparlett.jl")
include("aaa.jl")
include("ratkrylov.jl")

end # module
