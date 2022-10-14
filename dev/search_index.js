var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LinearEquations","category":"page"},{"location":"#LinearEquations","page":"Home","title":"LinearEquations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LinearEquations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LinearEquations]","category":"page"},{"location":"#LinearEquations.genLinearEquation-Tuple{Any, Any, Any, Any, Vararg{Any}}","page":"Home","title":"LinearEquations.genLinearEquation","text":"genLinearEquation(linearEquationName, sz, targetEquation, paramDisc, args...)\n\nget the time evolution form of an Master equation as linear equation form\n\nArguments:\n\nlinearEquationName: The name of linear equation you want to generate\nsz: The size of the system, number of energy levels\ntargetEquation: The target equation you want to linearized\nparamDisc: External code to convert the array p to parameters required\nargs: parameters required\n\nExamples:\n\njulia> function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)\n\tΓ::Float64 = 2π*5.746e6  # Decay Rate of Excited State\n    μ::Float64 = 2π*0.7e6\n    H = zeros(3, 3) * 0im\n    H[1, 1] = μ * B\n    H[2, 2] = -μ * B + δ\n    H[1, 3] = Ω / 2 * rabiRatio\n    H[2, 3] = Ω / 2\n    H[3, 1] = Ω' / 2 * rabiRatio\n    H[3, 2] = Ω' / 2\n    H[3, 3] = -1im/2 * Γ\n    return H\nend\njulia> function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, Ω::Float64, δ::Float64)\n\tΓ::Float64 = 2π*5.746e6  # Decay Rate of Excited State\n\trabiRatio::Float64 = 1\n\tΩc = Ω*(1+0im)\n    H = _getHamiltonian(B, Ωc, δ, rabiRatio)\n    ∂tρ = -1im * (H*ρ - ρ*H')\n    ∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)\n    ∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)\n    ∂tρ\nend\njulia> paramDisc = quote\n\t\t\tB0 = p[1]\n\t\t\tampB = p[2]\n\t\t\tphaseB = p[3]\n\t\t\tωB = p[4]\n\t\t\tprepare = p[5]\n\t\t\tB = (t > prepare) ? (ampB*sin(ωB*t+phaseB)+B0) : B0\n\t\t\tωΩ = p[6]\n\t\t\tτ = p[7]\n\t\t\tΩamp = p[8]\n\t\t\tT = 2π/ωΩ # 脉冲周期\n\t\t\tlowPulse = T - τ\n\t\t\tΩ = (t < prepare) | ((t>prepare) & ((t-prepare-τ/2)%T>lowPulse)) ? Ωamp : 0\n\t\t\tδ = p[9]\n\t\t\tend\njulia> expr = genLinearEquation(:master_equation!, 3, _∂tρ, paramDisc, :B, :Ω, :δ)\n\n\n\n\n\n","category":"method"},{"location":"#LinearEquations.getCoefficient-Tuple{Any, Any, Any}","page":"Home","title":"LinearEquations.getCoefficient","text":"getCoefficient(sz, argSize, masterFun)\n\nget the coefficient of an Master equation\n\nArguments:\n\ncoe: The coefficient dictionary.\nsz: The size of the system, number of energy levels\nmasterFun: The master equation\n\nExamples:\n\njulia> function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)\n\tΓ::Float64 = 2π*5.746e6  # Decay Rate of Excited State\n    μ::Float64 = 2π*0.7e6\n    H = zeros(3, 3) * 0im\n    H[1, 1] = μ * B\n    H[2, 2] = -μ * B + δ\n    H[1, 3] = Ω / 2 * rabiRatio\n    H[2, 3] = Ω / 2\n    H[3, 1] = Ω' / 2 * rabiRatio\n    H[3, 2] = Ω' / 2\n    H[3, 3] = -1im/2 * Γ\n    return H\nend\njulia> function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, Ω::Float64, δ::Float64)\n\tΓ::Float64 = 2π*5.746e6  # Decay Rate of Excited State\n\trabiRatio::Float64 = 1\n\tΩc = Ω*(1+0im)\n\tH = _getHamiltonian(B, Ωc, δ, rabiRatio)\n\t∂tρ = -1im * (H*ρ - ρ*H')\n\t∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)\n\t∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)\n\t∂tρ\nend\njulia> getCoefficient(3, 3, _∂tρ)\n\n\n\n\n\n","category":"method"},{"location":"#LinearEquations.vec2ρ-Tuple{Vector{Float64}}","page":"Home","title":"LinearEquations.vec2ρ","text":"vec2ρ(vec::Vector{Float64})::Matrix{ComplexF64}\n\nReconstruct density matrix from a vector, please refer ρ2vec\n\nArguments\n\nvec::Vector{Float64}: Vectorized of density matrix\n\nExamples\n\njulia> begin\n\tρ = zeros(3, 3) * 1im\n\tρ[1, 1] = 1/2\n\tρ[2, 2] = 1/2\n\tρ[1, 2] = 1 + 2im\n\tρ[2, 1] = 1 - 2im\n\tρ[3, 3] = 0\nend\njulia> vec = ρ2vec(ρ)\njulia> ρ2 = vec2ρ(vec)\njulia> ρ2 - ρ\n\n\n\n\n\n","category":"method"},{"location":"#LinearEquations.ρ2vec-Tuple{Matrix{ComplexF64}}","page":"Home","title":"LinearEquations.ρ2vec","text":"ρ2vec(ρ::Matrix{ComplexF64})::Vector{Float64}\n\nVectorize the density matrix ho, please refer vec2ρ\n\nArguments\n\nρ::Matrix{ComplexF64}: Density matrix\n\nExamples:\n\njulia> begin\n\tρ = zeros(3, 3) * 1im\n\tρ[1, 1] = 1/2\n\tρ[2, 2] = 1/2\n\tρ[1, 2] = 1 + 2im\n\tρ[2, 1] = 1 - 2im\n\tρ[3, 3] = 0\nend\njulia> vec = ρ2vec(ρ)\n\n\n\n\n\n","category":"method"}]
}
