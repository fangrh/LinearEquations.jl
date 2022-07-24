using LinearEquations
using Test

@testset "LinearEquations.jl" begin
    # Write your tests here.
    # function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)
    #     Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
    #     μ::Float64 = 2π*0.7e6
    #     H = zeros(3, 3) * 0im
    #     H[1, 1] = μ * B
    #     H[2, 2] = -μ * B + δ
    #     H[1, 3] = Ω / 2 * rabiRatio
    #     H[2, 3] = Ω / 2
    #     H[3, 1] = Ω' / 2 * rabiRatio
    #     H[3, 2] = Ω' / 2
    #     H[3, 3] = -1im/2 * Γ
    #     return H
    # end
    
    # function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, Ω::Float64, δ::Float64)
    #     Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
    #     rabiRatio::Float64 = 1
    #     Ωc = Ω*(1+0im)
    #     H = _getHamiltonian(B, Ωc, δ, rabiRatio)
    #     ∂tρ = -1im * (H*ρ - ρ*H')
    #     ∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)
    #     ∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)
    #     ∂tρ
    # end
    # @linearequation master_equation! 3 _∂tρ paramDisc B Ω δ
    function foo(args...)
        println("Hello World!")
        a = Vector([args...])
        println(a)
    end
    @testmacro foo
end
