module DynamicProgramming

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
    i    = 1
    Θᵢ   = keys(V₀)
    Vᵢ   = V₀
    Vᵢ₋₁ = Dict([(θ, -Inf) for θ ∈ Θᵢ])

    # Iterate
    while (maximum(abs.(get.(Ref(Vᵢ), Θᵢ, missing) - get.(Ref(Vᵢ₋₁), Θᵢ, missing))) >= tolerance && i <= maxiter)
        # Save the state
        Vᵢ₋₁ = Vᵢ

        # Do the grid search
        for θ ∈ Θᵢ
            # Initialize search value
            Vᵢ[θ] = -Inf

            # Loop over potential choices
            for ϑ ∈ Θᵢ
                # Evaluate choice value
                Vₖₗ = u(θ, ϑ) + β*Vᵢ₋₁[ϑ]
                Vᵢ[θ] = (Vₖₗ > Vᵢ[θ]) ? Vₖₗ : Vᵢ[θ]
            end
        end
    end

    θ₊ = Dict([(θ, θ) for θ ∈ Θᵢ])
    V₊ = Vᵢ
    # Extract policy function
    for θ ∈ Θᵢ
        # Initialize search value
        V₊[θ] = -Inf

        # Loop over potential choices
        for ϑ ∈ Θᵢ
            # Evaluate choice value
            Vₖₗ = u(θ, ϑ) + β*Vᵢ[ϑ]
            θ₊[θ] = (Vₖₗ > V₊[θ]) ? ϑ : θ₊[θ]
        end
    end

    return Vᵢ, θ₊
end

u(k, k̃, a) = log(k^a - k̃)
uₐ(k, k̃) = u(k, k̃, 0.1)
β = 0.9
V₀ = Dict(0.1 => 0.1, 0.2 => 0.2, 0.3 => 0.3)

solveVFI(uₐ, β, V₀)

end