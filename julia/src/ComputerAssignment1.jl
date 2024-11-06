module ComputerAssignment1

using Revise
using Plots
using LinearAlgebra

#== Problem 1 ==#
includet("DynamicProgramming/DynamicProgramming.jl")
import .DynamicProgramming

# Parameters
α = 1.5
β = 0.9
δ = 0.9

# Exogenous
z = 1

# Functions
f(k) = k^α # Production function
u(c) = c > 0 ? log(c) : -Inf # Instantaneous utility
c(k, κ) = f(k) + k*(1 - δ) - κ # Consumption policy function
U(k, κ) = u(c(k, κ)) # Derived utility

# Initialization
k̲ = 0.81
k̄ = 2.00
K  = k̲:0.001:k̄
V₀ = Dict([(k, U(k, k)/(1-β)) for k ∈ K])

# Solve the problem
V̂, k̂ = DynamicProgramming.solveVFI(U, β, V₀; maxiter=100)

# Plot value function
k₀ = k̲:0.1:k̄
T  = 15
kₜ, Vₜ = DynamicProgramming.policyTrace(T, k̂, k₀, V̂)
p = plot()
for k ∈ k₀
    plot!(p, 1:T, kₜ, palette=:Dark2_5, legend=false)
end
display(p)


# Problem 5
include("Roots/Roots.jl")
import .Roots


end