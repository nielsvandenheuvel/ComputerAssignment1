module DynamicProgramming

using Printf
using LinearAlgebra

# Type alias for working input
const FunctionOnGrid = AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real}

"""
# Neoclassical Growth Model Solver Using Value Function Iteration Algorithm
Method to find value- and policy function.

## Positional Arguments:
- __u__ : : Function | Transition utility function; that is, utility from going from one state to another, 
with the arguments positional and in this order.
- __β__ : : Real | Discount factor.
- __V₀__ : : FunctionOnGrid | The initial guess for the value function.

## Named Arguments
- __monotone__ : : Bool | Exploit monotonicity of the value function to refine the grid.
- __concave__ : : Bool | Exploit concavity of the value function to refine the grid.
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
- FunctionOnGrid | Value function.
- FunctionOnGrid | Policy function.

_P.S. FunctionOnGrid is a type alias for AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real}._
"""
function solveVFI(
    u::Function,
    β::Real,
    V₀::FunctionOnGrid;
    monotone::Bool = false,
    concave::Bool = false,
    tolerance::Real = 1e-5,
    maxiter::Real = 1e5
)
    @assert β >= 0 && β <= 1 "The discount factor must be in the interval [0,1]."

    # Initialize
    i   = 1
    Vᵢ  = V₀
    πᵢ  = Dict([(θ, θ) for θ ∈ eachindex(Vᵢ)])
    ΔVᵢ = tolerance

    # Iterate
    while (ΔVᵢ >= tolerance && i <= maxiter)
        # Measure performance
        start_time = time_ns()

        # Save the state
        ΔVᵢ  = 0
        Vᵢ₋₁ = copy(Vᵢ)
        
        # Do the grid search
        for θ ∈ keys(Vᵢ)
            # Loop over potential choices
            for ϑ ∈ keys(Vᵢ)
                # No need to check if new value is -Inf
                if (Vᵢ₋₁[ϑ] > -Inf)
                    # Evaluate choice value
                    V̂ = u(θ, ϑ) + β*Vᵢ₋₁[ϑ]
                    if (V̂ > Vᵢ[θ])
                        Vᵢ[θ] = V̂
                        πᵢ[θ] = ϑ
                    end
                end
            end

            # Update slack
            if (ΔVᵢ < Vᵢ[θ] - Vᵢ₋₁[θ])
                ΔVᵢ = Vᵢ[θ] - Vᵢ₋₁[θ]
            end
        end

        # Update
        i += 1
        end_time = time_ns()
        if (i % (maxiter/10) == 0)
            println("Progress: $(round(i/maxiter*100))% | Slack: $(@sprintf("%.2e", ΔVᵢ)) | ETA: $(round((end_time-start_time)*(maxiter-i)/1e9/60, digits=2))min")
        end
    end

    return Vᵢ, πᵢ
end

"""
# Tracer for Policy Function
Method to trace out a policy function for a number of steps. _TD: inverse tracing does not work yet._.

## Positional Arguments:
- __T__ : : Integer | Number of steps to trace.
- __π__ : : FunctionOnGrid | Policy function.
- __θ₀__ : : Real | The initial state.

## Returns:
## Returns:
- Vector | Policy trace.

_P.S. FunctionOnGrid is a type alias for AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real}._
"""
function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::Real) :: Vector
    if (T == 0)
        return θ₀
    else
        # Initialize path vector
        Π = zeros(abs(T))
        
        # Allow to go back in time
        if (T > 0)
            # Intialize
            Π[1] = π[θ₀]

            # Trace out path
            for t ∈ 2:T
                Π[t] = π[Π[t-1]]
            end
        else
            # Intialize
            Π[-T] = [θ for θ ∈ keys(π) if θ .== θ₀][1]

            # Trace out path
            for t ∈ (T+1):-1
                Π[-t] = [θ for θ ∈ keys(π) if θ .== Π[1-t]][1]
            end
        end
    end

    return Π
end

"""
# Tracer for Policy Function
Wrapper to trace out value corresponding to policy trace.

## Additonal Positional Arguments:
- __v__ : : FunctionOnGrid | Value function.

_P.S. FunctionOnGrid is a type alias for AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real}._
"""
function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::Real, v::FunctionOnGrid) :: Tuple{Vector, Vector}
    # Obtain policy trace
    Π = policyTrace(T, π, θ₀)

    # Initialize path vector
    V = zeros(abs(T))

    # Trace out path
    for t ∈ 1:T
        V[t] = v[Π[t]]
    end

    return Π, V
end

"""
## Tracer for Policy Function
Wrapper for tracing multiple initial states.
"""
function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::AbstractVector) :: Matrix
    # Initialize path matrix
    Π = zeros(abs(T), length(θ₀))

    # Fill
    for (k, θ) ∈ enumerate(θ₀)
        Π[1:T, k] = policyTrace(T, π, θ)
    end

    return Π
end

"""
## Tracer for Policy Function
Wrapper for tracing multiple initial states.
"""
function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::AbstractVector, v::FunctionOnGrid) :: Tuple{Matrix, Matrix}
    # Initialize path matrix
    Π = zeros(abs(T), length(θ₀))
    V = zeros(abs(T), length(θ₀))

    # Fill
    for (k, θ) ∈ enumerate(θ₀)
        Π[1:T, k], V[1:T, k] = policyTrace(T, π, θ, v)
    end

    return Π, V
end

end