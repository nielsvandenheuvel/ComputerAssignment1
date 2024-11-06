module ComputerAssignment1

using Revise
using Plots

#== Problem 1 ==#
includet("DynamicProgramming/DynamicProgramming.jl")
import .DynamicProgramming

# Parameters
α = 0.5
β = 0.9
δ = 0.9
T = 100

# Exogenous
z = 1

# Functions
F(k) = k^α # Production function
f(k) = α*k^(α-1) # Marginal production
u(c) = c > 0 ? log(c) : -Inf # Instantaneous utility
c(k, κ) = k*(1 - δ + z*f(k)) - κ > 0 ? k*(1 - δ + z*f(k)) - κ : 1E-5 # Consumption policy function
ũ(k, κ) = u(c(k, κ)) # Derived utility

# Initialization
K  = 0.01:0.01:01
V₀ = Dict([(k, ũ(k, k)/(1-β)) for k ∈ K])

# Solve the problem
V̂, k̂ = DynamicProgramming.solveVFI(ũ, β, V₀; maxiter=100)

# Plot value function
k₀ = 0.7
T = 100
k̂₀ = zeros(T,1)
k̂₀[1] = k₀
V̂₀ = zeros(T,1)
V̂₀[1] = V̂[k₀]
for t ∈ 2:T
    k̂₀[t] = k̂[k̂₀[t-1]]
    V̂₀[t] = V̂[k̂₀[t]]
end
plot(1:T, k̂₀)
plot(1:T, V̂₀)


# Problem 5
include("Roots/Roots.jl")
import .Roots


end