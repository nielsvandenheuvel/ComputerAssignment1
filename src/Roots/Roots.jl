module Roots

using LinearAlgebra

include("bisection.jl")
include("newton.jl")
include("quasi_newton.jl")
include("inverse_quadratic_interpolation.jl")
include("contraction.jl")

end