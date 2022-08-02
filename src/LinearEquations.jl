module LinearEquations
export genLinearEquation, vec2ρ, ρ2vec
# export @testmacro

using Printf


"""
@testmacro funName

This macro is used for test some features
"""
# macro testmacro(funName)
# 	args = [1, 2, 3]
# 	# expr1 = quote
# 	# 	$funName($args...)
# 	# end
# 	# esc(expr1)
# 	b = esc(:($funName($args...)))
# 	println(b)
# end

"""
ρ2vec(ρ::Matrix{ComplexF64})::Vector{Float64}

Vectorize the density matrix `\rho`, please refer [vec2ρ](@Ref)

# Arguments
- `ρ::Matrix{ComplexF64}`: Density matrix

# Examples:

```julia-repl
julia> begin
	ρ = zeros(3, 3) * 1im
	ρ[1, 1] = 1/2
	ρ[2, 2] = 1/2
	ρ[1, 2] = 1 + 2im
	ρ[2, 1] = 1 - 2im
	ρ[3, 3] = 0
end
julia> vec = ρ2vec(ρ)
```
"""
function ρ2vec(ρ::Matrix{ComplexF64})::Vector{Float64}
	sz = size(ρ)[1]
	vec = zeros(sz*sz)
	for i = 1:sz
		vec[i] = real(ρ[i, i])
	end
	index = sz + 1
	for i = 1:(sz-1)
		for j = (i+1) : sz
			vec[index] = real(ρ[i, j])
			vec[index+1] = imag(ρ[i, j])
			index += 2
		end
	end
	vec
end;

"""
vec2ρ(vec::Vector{Float64})::Matrix{ComplexF64}

Reconstruct density matrix from a vector, please refer [ρ2vec](@Ref)

# Arguments

- `vec::Vector{Float64}`: Vectorized of density matrix

# Examples
```julia-repl
julia> begin
	ρ = zeros(3, 3) * 1im
	ρ[1, 1] = 1/2
	ρ[2, 2] = 1/2
	ρ[1, 2] = 1 + 2im
	ρ[2, 1] = 1 - 2im
	ρ[3, 3] = 0
end
julia> vec = ρ2vec(ρ)
julia> ρ2 = vec2ρ(vec)
julia> ρ2 - ρ
```
"""
function vec2ρ(vec::Vector{Float64})::Matrix{ComplexF64}
	sz = floor(Int64, sqrt(size(vec)[1]))
	ρ = zeros(sz, sz) * 1im
	for i = 1:sz
		ρ[i, i] = vec[i]
	end
	index = sz + 1
	for i = 1:(sz - 1)
		for j = (i+1) : sz
			ρ[i, j] = vec[index] + 1im*vec[index + 1]
			ρ[j, i] = vec[index] - 1im*vec[index + 1]
			index += 2
		end
	end
	ρ
end;


"""
genLinearEquation(linearEquationName, sz, targetEquation, paramDisc, args...)

get the time evolution form of an Master equation as linear equation form

# Arguments:

- `linearEquationName`: The name of linear equation you want to generate
- `sz`: The size of the system, number of energy levels
- `targetEquation`: The target equation you want to linearized
- `paramDisc`: External code to convert the array p to parameters required
- `args`: parameters required
# Examples:

```julia-repl
julia> function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)
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
julia> function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, Ω::Float64, δ::Float64)
	Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
	rabiRatio::Float64 = 1
	Ωc = Ω*(1+0im)
    H = _getHamiltonian(B, Ωc, δ, rabiRatio)
    ∂tρ = -1im * (H*ρ - ρ*H')
    ∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)
    ∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)
    ∂tρ
end
julia> paramDisc = quote
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
julia> expr = genLinearEquation(:master_equation!, 3, _∂tρ, paramDisc, :B, :Ω, :δ)
```
"""
function genLinearEquation(linearEquationName, sz, targetEquation, paramDisc, args...)
	sysSize = sz * sz
	argSize = length(args)
	exprs = []
	dρExprStr = [@sprintf "du[%d] = " i for i in 1:sysSize] # 演化的结果
	for m = 1:argSize
		argVec = zeros(argSize)
		argVec[m] = 1
		for j = 1:sysSize
			vec = zeros(sysSize)
			vec[j] = 1
			ρ = vec2ρ(vec)
			dρ = targetEquation(ρ, argVec...)
			dρ0 = targetEquation(ρ, zeros(argSize)...)
			dρVec = ρ2vec(dρ-dρ0)	
			for k = 1:sysSize
				if dρVec[k] > 0 # 如果不等于零，那么主方程就要加上相应的值
					# println(args[m])
					if dρExprStr[k] != (@sprintf "du[%d] = " k)
						addStr = @sprintf "+ %f * %s * u[%d]" dρVec[k] args[m] j
					else
						addStr = @sprintf " %f * %s * u[%d]" dρVec[k] args[m] j
					end
					dρExprStr[k] = dρExprStr[k] * addStr
				elseif dρVec[k] < 0
					addStr = @sprintf " %f * %s * u[%d]" dρVec[k] args[m] j
					dρExprStr[k] = dρExprStr[k] * addStr
				end
			end
		end
	end
	for j = 1:sysSize
		vec = zeros(sysSize)
		vec[j] = 1
		ρ = vec2ρ(vec)
		argVec = zeros(argSize)
		dρ = targetEquation(ρ, argVec...)
		dρVec = ρ2vec(dρ)
		# println(dρVec)
		# println(dρ)
		# println(vec)
		# println(argVec)
		
		for k = 1:sysSize
			if dρVec[k] > 0 # 如果不等于零，那么主方程就要加上相应的值
				if dρExprStr[k] != (@sprintf "du[%d] = " k)
					addStr = @sprintf "+ %f * u[%d]" dρVec[k] j
				else
					addStr = @sprintf " %f * u[%d]" dρVec[k] j
				end
				dρExprStr[k] = dρExprStr[k] * addStr
			elseif dρVec[k] < 0
				addStr = @sprintf " %f * u[%d]" dρVec[k] j
				dρExprStr[k] = dρExprStr[k] * addStr
			end
		end
	end
	# for i=1:sysSize
	# 	println(dρExprStr[i])
	# end
	push!(exprs, paramDisc)
	for l = 1:sysSize
		if dρExprStr[l] != (@sprintf "du[%d] = " l)  # 如果最终没有加入，那么就不写进来
			expri = Meta.parse(dρExprStr[l])
			push!(exprs, expri)
		end
	end

	quote
		function $(linearEquationName)(du, u, p, t)
			$(exprs...)
		end
	end
end

"""
getCoefficient(sz, argSize, masterFun)

get the coefficient of an Master equation

# Arguments:

- `coe`: The coefficient dictionary.
- `sz`: The size of the system, number of energy levels
- `masterFun`: The master equation

# Examples:

```julia-repl
julia> function _getHamiltonian(B::Float64, Ω::ComplexF64, δ::Float64, rabiRatio::Float64)
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
julia> function _∂tρ(ρ::Matrix{ComplexF64}, B::Float64, ReΩ::Float64, ImΩ::Float64, δ::Float64)
	Γ::Float64 = 2π*5.746e6  # Decay Rate of Excited State
	Ω = ReΩ + 1im*ImΩ
	rabiRatio::Float64 = 1
    H = _getHamiltonian(B, Ω, δ, rabiRatio)
    ∂tρ = -1im * (H*ρ - ρ*H')
    ∂tρ[1, 1] = ∂tρ[1, 1] + ρ[3, 3] * Γ * rabiRatio^2 / (rabiRatio^2 + 1)
    ∂tρ[2, 2] = ∂tρ[2, 2] + ρ[3, 3] * Γ / (rabiRatio^2 + 1)
    ∂tρ
end
julia> getcoefficient(3, 4, _∂tρ)
```
"""
function getCoefficient(sz, argSize, masterFun)
	sysSize = sz^2
	coeM_arr = []
	#-----------------------------
	p = zeros(argSize)
	coeM0 = zeros(argSize, argSize) 
	for i = 1:sysSize
		u = zeros(sysSize)
		u[i] = 1.0
		ρ = vec2ρ(u)
		dρ = masterFun(ρ, p...)
		coeM0[:, i] = ρ2vec(dρ)
	end
	push!(coeM_arr, coeM0)
	#-----------------------------
	for a = 1:argSize
		p = zeros(argSize)
		p[a] = 1.0
		coeM = zeros(argSize, argSize) 
		for i = 1:sysSize
			u = zeros(sysSize)
			u[i] = 1.0
			ρ = vec2ρ(u)
			dρ = masterFun(ρ, p...)
			coeM[:, i] = ρ2vec(dρ)
		end
		push!(coeM_arr, coeM-coeM0)
	end
	return coeM_arr
end

end
