using DifferentialEquations, ForwardDiff, LinearAlgebra

abstract type DynamicalSystem end

mutable struct DiscreteDynamicalSystem <: DynamicalSystem
	rule::Function
	u::Vector{Float64}
end

mutable struct ContinuousDynamicalSystem <: DynamicalSystem
	rule::Function     # f(u) (works for ForwardDiff)
	rule_de::Function  # f(du,u,p,t) (works for DifferentialEquations) 
	u::Vector{Float64} #
	p::Vector{Float64} # for DifferentialEquations
end

function lyapunovspectrum(ds::ContinuousDynamicalSystem, N::Float64, Δt::Float64, k::Int64)
	# t∈[0,T] & T = Δt*N
	D = length(ds.u)
	if k > D error("k > D not allowed.") end
	Xprob = ODEProblem(ds.rule_de, ds.u, (0.0, Δt*N), ds.p)
	Xint = init(Xprob)

	function Y_evolution(du, u, p, t)
		Y = reshape(u, (D,k))
		J = ForwardDiff.jacobian(ds.rule, Xint.u)
		Ẏ = vec(J * Y)
		for i=1:length(du)
			du[i] = Ẏ[i]
		end
	end

	Yprob = ODEProblem(Y_evolution, vec(Matrix{Float64}(I, D, k)), (0.0, Δt*N))
	Yint = init(Yprob)

	λ = zeros(k)
	for _ in 1:N
		step!(Xint, Δt, true)
		step!(Yint, Δt, true)
		Y = reshape(Yint.u, (D,k))
		Q, R = qr(Y)
		Q = Matrix(Q)
		λ .+= log.(abs.(diag(R)))
		reinit!(Yint, vec(Q))
	end
	λ ./= (Δt*N)
	return λ
end


function sprott_rule(u)
	x, y, z = u
	return [
		y + 2*x*y + x*z,
		1 - 2*x^2 + y*z,
		x - x^2 - y^2
	]
end

function sprott_rule_de!(du, u, p, t)
	x, y, z = u
	du[1] = y + 2*x*y + x*z
	du[2] = 1 - 2*x^2 + y*z
	du[3] = x - x^2 - y^2
end

sprott_con = ContinuousDynamicalSystem(sprott_rule, sprott_rule_de!, [1.0,0.0,0.0], [])
sprott_dis = ContinuousDynamicalSystem(sprott_rule, sprott_rule_de!, [2.0,0.0,0.0], [])

ls_con = lyapunovspectrum(sprott_con, 1_000.0, 1.0, 3)
ls_dis = lyapunovspectrum(sprott_dis, 1_000.0, 1.0, 3)

println("Conservative: sum(λ) = $(round(sum(ls_con), digits=4)) ≈ 0\nDissipative:  sum(λ) = $(round(sum(ls_dis), digits=4)) < 0 and λ1 = $(round(ls_dis[1], digits=4)) > 0")
#=
Conservative: sum(λ) = -0.0009 ≈ 0
Dissipative:  sum(λ) = -0.1602 < 0 and λ1 = 0.0735 > 0
=#
