include("../src/Roots/Roots.jl")
import .Roots

@testset "Roots" begin
    @testset "Unidimensional" begin
        # Functions
        f = (x::Real) -> x      # Linear
        g = (x::Real) -> x^2    # Quadratic
        h = (x::Real) -> sin(x) # Sinusoid
        
        # Tests
        @testset "Bisection" begin
            @test Roots.solve(f; method="Bisection", a=-1, b=2) == 0
        end
    end
end