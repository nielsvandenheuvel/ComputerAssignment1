module DynamicProgramming

using Printf
using LinearAlgebra

# Type alias for working input
const FunctionOnGrid = AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real}

function solveVFI(
    u::Function,
    β::Real,
    V₀::FunctionOnGrid;
    monotone::Bool = false,
    concave::Bool = false,
    tolerance::Real = 10E-5,
    maxiter::Real = 10E5
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

function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::Real; plot::Bool = false)
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

    # Plot results
    if (plot)

    end

    return Π
end

function policyTrace(T::Integer, π::FunctionOnGrid, θ₀::AbstractVector{<:Real})
    # Initialize path matrix
    Π = zeros(abs(T), length(θ₀))

    # Fill
    for (k, θ) ∈ enumerate(θ₀)
        Π[1:T, k] = policyTrace(T, π, θ)
    end

    return Π
end

function policyTrace

end