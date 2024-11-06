"""
# Contraction Root Finding Algorithm
Method to find roots using the contraction property.

## Positional Arguments:
- __f__ : : Function | The function for which you want to find the roots.
- __x₀__ : : Real | The initial guess.

## Named Arguments
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
The found root.
"""
function solveContraction(f::Function, x₀::Union{<:Real, <:AbstractArray{<:Real}}; tolerance::Real = 10E-5, maxiter::Real = 10E5)
    @assert tolerance > 0 "The convergence criterion should be positive."

    # Initialize
    i  = 1
    xᵢ = x₀
    fᵢ = f(x₀)
    
    # Iterate
    while (norm(fᵢ) >= tolerance && i <= maxiter && abs(xᵢ-f(xᵢ)) < 1/tolerance)
        # Update
        fᵢ = f(xᵢ)
        xᵢ = fᵢ
        i += 1
    end

    return xᵢ
end