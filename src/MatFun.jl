__precompile__()

module MatFun

using Base.LinAlg, TaylorSeries

export schurparlett, ratkrylov, ratkrylovf

include("schurparlett.jl")
include("ratkrylov.jl")

end # module
