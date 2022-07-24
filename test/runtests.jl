using LinearEquations
using Test

@testset "LinearEquations.jl" begin
    # Write your tests here.
    function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)
        Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
        μ::Float64 = 2π*0.7e6
        H = zeros(3, 3) * 0im
        H[1, 1] = μ * B
        H[2, 2] = -μ * B + δ
        H[1, 3] = Ω / 2 * rabiRatio
        H[2, 3] = Ω / 2
        H[3, 1] = Ω' / 2 * rabiRatio
        H[3, 2] = Ω' / 2
        H[3, 3] = -1im/2 * Γ
        return H
    end
    
    function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, Ω::Float64, δ::Float64)
        Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
        rabiRatio::Float64 = 1
        Ωc = Ω*(1+0im)
        H = _getHamiltonian(B, Ωc, δ, rabiRatio)
        ∂tρ = -1im * (H*ρ - ρ*H')
        ∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)
        ∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)
        ∂tρ
    end
    paramDisc = quote
        B0 = p[1]
        ampB = p[2]
        phaseB = p[3]
        ωB = p[4]
        prepare = p[5]
        B = (t > prepare) ? (ampB*sin(ωB*t+phaseB)+B0) : B0
        ωΩ = p[6]
        τ = p[7]
        Ωamp = p[8]
        T = 2π/ωΩ # 脉冲周期
        lowPulse = T - τ
        Ω = (t < prepare) | ((t>prepare) & ((t-prepare-τ/2)%T>lowPulse)) ? Ωamp : 0
        δ = p[9]
    end
    # @linearequation master_equation! 3 _∂tρ paramDisc B Ω δ
    expr1 = genLinearEquation(:master_equation!, 3, _∂tρ, paramDisc, :B, :Ω, :δ)
    println(expr1)
    eval(expr1)
    # @linearequation2 master_equation! 3 _∂tρ paramDisc B Ω δ

    # function foo(args...)
    #     println("Hello World!")
    #     a = Vector([args...])
    #     println(a)
    #     b = a[1] + a[2] + a[3]
    #     return b
    # end
    # @testmacro eval(foo)
end
