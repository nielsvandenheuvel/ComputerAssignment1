"""
# Secant Root Finding Algorithm
Method to find roots using the secant algorithm.

## Positional Arguments:
- __f__ : : Function | The function for which you want to find the roots.
- __x₁__ : : Real | The second initial value and also first guess.
- __x₀__ : : Real | The first initial value.

## Named Arguments:
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
The found root.
"""
function solveSecant(f::Function,
    x₁::Real,
    x₀::Real;
    tolerance::Real = 10E-5,
    maxiter::Real = 10E5
) :: Real
    # Initialize
    i    = 1
    xᵢ   = x₁
    xᵢ₋₁ = x₀

    # Iterate
    while (abs(f(xᵢ)) >= tolerance && i <= maxiter)
        # Save the state
        sᵢ = xᵢ₋₁
        # Update
        xᵢ  -= (xᵢ - xᵢ₋₁)/(f(xᵢ) - f(xᵢ₋₁))*f(xᵢ)
        xᵢ₋₁ = sᵢ
        i   += 1
    end

    return xᵢ
end

"""
# Broyden's Root Finding Algorithm
Method to find roots using Broyden's algorithm with the Sherman-Morrison formula for the approximate Jacobian.

## Positional Arguments:
- __f__ : : Function | The function for which you want to find the roots.
- __x₀__ : : Vector{Real} | The initial value.
- __S₁__ : : Matrix{Real} | The initial approximate inverse Jacobian.

## Named Arguments:
- __scale__ : : Real | The scaling factor to control the step size.
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
The found root.
"""
function solveBroyden(
    f::Function,
    x₀::AbstractVector{<:Real},
    S₁::Union{AbstractMatrix{<:Real}, Nothing} = nothing;
    scale::Real = 1,
    tolerance::Real = 10E-5,
    maxiter::Real = 10E5
) :: Vector{<:Real}

    if (S₁ !== nothing)
        @assert size(S₁, 1) == length(x₀) "The inverse Jacobian should have the same number of rows as the initial vectors."
        @assert size(S₁, 1) == size(S₁, 2) "The inverse Jacobian should be square."
    else
        S₁ = Matrix{Real}(I, length(x₀))
    end

    # Initialize
    i  = 1
    xᵢ = xᵢ₋₁ = x₀
    Sᵢ = Sᵢ₋₁ = S₁

    # Iterate
    while (abs(f(xᵢ)) >= tolerance && i <= maxiter)
        # Save the state
        sᵢ  = xᵢ₋₁
        # Update
        xᵢ  -= scale*Sᵢ*f(xᵢ)
        Δxᵢ  = xᵢ - xᵢ₋₁
        Δfᵢ  = f(xᵢ) - f(xᵢ₋₁)
        Sᵢ  += (Δxᵢ - Sᵢ₋₁*Δfᵢ)/(norm(Δfᵢ)^2)*transpose(Δfᵢ)
        xᵢ₋₁ = sᵢ
        i   += 1
    end

    return vec(xᵢ)
end