"""
# Univariate Inverse Quadratic Interpolation Root Finding Algorithm
Method to find roots using the inverse quadratic interpolation algorithm.

## Positional Arguments:
- __f__ : : Function | The function for which you want to find the roots.
- __x₂__ : : Real | The third initial value and also first guess.
- __x₁__ : : Real | The second initial value.
- __x₀__ : : Real | The first initial value.

## Named Arguments
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
The found root.
"""
function solveIQI(f::Function, x₂::Real, x₁::Real, x₀::Real; tolerance::Real = 10E-5, maxiter::Real = 10E5)
    @assert tolerance > 0 "The convergence criterion should be positive."

    # Initialize
    i    = 1
    xᵢ   = x₂
    xᵢ₋₁ = x₁
    xᵢ₋₂ = x₀

    # Iterate
    while (abs(f(xᵢ)) >= tolerance && i <= maxiter)
        # Save the state
        sᵢ = (xᵢ₋₁, xᵢ₋₂)
        # Update
        xᵢ   = (
            (f(xᵢ₋₁)f(xᵢ))/((f(xᵢ₋₂)-f(xᵢ₋₁))*(f(xᵢ₋₂)-f(xᵢ)))xᵢ₋₂ +
            (f(xᵢ₋₂)f(xᵢ))/((f(xᵢ₋₁)-f(xᵢ₋₂))*(f(xᵢ₋₁)-f(xᵢ)))xᵢ₋₁ +
            (f(xᵢ₋₁)f(xᵢ₋₂))/((f(xᵢ)-f(xᵢ₋₂))*(f(xᵢ)-f(xᵢ₋₁)))xᵢ
        )
        xᵢ₋₁ = sᵢ[1]
        xᵢ₋₂ = sᵢ[2]
        i   += 1
    end

    return xᵢ
end