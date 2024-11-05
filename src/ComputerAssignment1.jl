using Plots

# Problem 5
include("Roots/Roots.jl")
import .Roots

# The function
f(x) = sin(2π*x) - (1/2)x

# Solve for different intial values
rootsBisection(δ)    = Roots.solveBisection(f, -δ, δ)
rootsSecant(x₀)      = Roots.solveSecant(f, x₀, x₀+0.3)
rootsContraction(x₀) = Roots.solveContraction(f, x₀, maxiter=100000)

plot(-2:0.01:2, rootsContraction, seriestype=:scatter)