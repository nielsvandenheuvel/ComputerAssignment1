module DynamicProgramming

using Printf

function solveVFI(
    u::Function,
    β::Real,
    V₀::AbstractDict{<:Union{<:Real, <:AbstractVector{<:Real}}, <:Real};
    monotone::Bool = false,
    concave::Bool = false,
    tolerance::Real = 10E-5,
    maxiter::Real = 10E5
)
    @assert β >= 0 && β <= 1 "The discount factor must be in the interval [0,1]."

    # Initialize
    i   = 1
    Vᵢ  = V₀
    θᵢ  = Dict([(θ, θ) for θ ∈ eachindex(Vᵢ)])
    ΔVᵢ = tolerance

    # Iterate
    while (ΔVᵢ >= tolerance && i <= maxiter)
        # Save the state
        Vᵢ₋₁ = Vᵢ

        # Do the grid search
        for θ ∈ keys(Vᵢ)
            # Initialize search value
            Vᵢ[θ] = -Inf
            # Loop over potential choices
            for ϑ ∈ keys(Vᵢ)
                # Evaluate choice value
                Vₖₗ = u(θ, ϑ) + β*Vᵢ₋₁[ϑ]
                if (Vₖₗ > Vᵢ[θ])
                    Vᵢ[θ] = Vₖₗ
                    θᵢ[θ] = ϑ
                end
            end

            # Update slack
            if (ΔVᵢ < Vᵢ[θ] - Vᵢ₋₁[θ])
                ΔVᵢ = Vᵢ[θ] - Vᵢ₋₁[θ]
            end
        end

        # Update
        i += 1
        if (i % (maxiter/10) == 0)
            println("Iteration $i/$maxiter: slack is $(@sprintf("%.2e", ΔVᵢ))")
        end
    end

    return Vᵢ, θᵢ
end

end