module LinearEquations
export @linearequation

using Printf

"""
@targetEquation(targetEquationName, sz, targetEquation, paramDisc, args...)

get the time evolution form of an Master equation

# Arguments:

- `targetEquationName`: The name of linear equation you want to generate
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
julia> @targetEquation master_equation! 3 _∂tρ paramDisc B Ω δ
```
"""
macro linearequation(targetEquationName, sz, targetEquation, paramDisc, args...)
	sysSize = eval(sz) * eval(sz)
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
			dρ = eval(targetEquation)(ρ, argVec...)
			dρ0 = eval(targetEquation)(ρ, zeros(argSize)...)
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
		dρ = eval(targetEquation)(ρ, argVec...)
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
	push!(exprs, eval(paramDisc))
	for l = 1:sysSize
		if dρExprStr[l] != (@sprintf "du[%d] = " l)  # 如果最终没有加入，那么就不写进来
			expri = Meta.parse(dρExprStr[l])
			push!(exprs, expri)
		end
	end

	quote
		function $(esc(targetEquationName))(du, u, p, t)
			$(exprs...)
		end
	end
end

end
