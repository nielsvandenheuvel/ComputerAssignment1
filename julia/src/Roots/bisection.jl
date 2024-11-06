"""
# Bissection Root Finding Algorithm
Method to find roots using the bissection algorithm.

## Positional Arguments:
- __f__ : : Function | The function for which you want to find the roots.
- __a__ : : Real | The initial lower bound.
- __b__ : : Real | The initial upper bound.

## Named Arguments
- __tolerance__ : : Real | Stopping tolerance.
- __maxiter__ : : Real | Maximum number of iterations.

## Returns:
The found root.
"""
function solveBisection(f::Function, a::Real, b::Real; tolerance::Real = 10E-5, maxiter::Real = 10E5)
    @assert a < b "The lower bound should be strictly smaller than the upper bound."
    @assert f(a) >= 0 >= f(b) || f(a) <= 0 <= f(b) "f(a) and f(b) need to of opposite sign."

    # Check if boundaries are solutions
    if (abs(f(a)) < tolerance)
        return a
    end
    if (abs(f(b)) < tolerance)
        return b
    end

    # Initialize
    i  = 1
    aᵢ = a
    bᵢ = b
    xᵢ = (aᵢ + bᵢ)/2
    
    # Iterate
    while (abs(f(xᵢ)) >= tolerance && i <= maxiter)
        # Reset the interval
        if (sign(f(xᵢ)) == sign(f(aᵢ)))
            aᵢ = xᵢ
        else
            bᵢ = xᵢ
        end

        # Update
        xᵢ = (aᵢ + bᵢ)/2
        i += 1
    end

    return xᵢ
end