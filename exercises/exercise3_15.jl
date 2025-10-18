using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots
using PredefinedDynamicalSystems

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

function lyapunovspectrum_convergance(ds::ContinuousDynamicalSystem, N::Float64, Δt::Float64, k::Int64)
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

	λN = Vector{Vector{Float64}}()
	λ = zeros(k)
	for i in 1:N
		step!(Xint, Δt, true)
		step!(Yint, Δt, true)
		Y = reshape(Yint.u, (D,k))
		Q, R = qr(Y)
		Q = Matrix(Q)
		λ .+= log.(abs.(diag(R)))
		reinit!(Yint, vec(Q))
		push!(λN, λ./(Δt*i))
	end
	return λN
end

function lyapunovspectrum_convergance_time(lyapunovspectrum_convergance::Vector{Vector{Float64}}, Δt; ϵ=0.1)::Float64
	δ = [ls .- lyapunovspectrum_convergance[end] for ls in lyapunovspectrum_convergance]
	return Δt * findfirst(x->all(y->y<=ϵ, abs.(x)), δ)
end


let u0s = PredefinedDynamicalSystems.henonheiles_ics(0.13, 15), T = 1_250.0, Δt = 1.0, pc_ϵ=0.1
	poincare_section = Vector{Tuple{Vector{Float64}, Float64}}()

	function henon_heiles_rule(u)
		x, y, vx, vy = u
		return [
			vx
			vy
			-x - 2*x*y
			-y - (x^2 - y^2)
		]
	end
	function henon_heiles_rule_de!(du, u, p, t)
		x, y, vx, vy = u
		du[1] = vx
		du[2] = vy
		du[3] = -x - 2*x*y
		du[4] = -y - (x^2 - y^2)
	end

	for (i, u0) in enumerate(u0s)
		println(vcat([" " for _ in 1:(length(string(length(u0s)))-length(string(i)))], ["$i", " of $(length(u0s))"])...)
		lo = ContinuousDynamicalSystem(henon_heiles_rule, henon_heiles_rule_de!, u0, [])
		lsN = lyapunovspectrum_convergance(lo, T, Δt, length(u0))
		lsct = lyapunovspectrum_convergance_time(lsN, Δt; ϵ=0.01)
		llsct = log(lsct)

		prob = ODEProblem(henon_heiles_rule_de!, u0, (0.0,T))
		sol  = solve(prob)
		
		fa = findall(x->(abs(x[1])<=pc_ϵ), sol.u)
		points = sol.u[fa]

		push!(poincare_section, [(p,llsct) for p in points]...)
	end

	scatter(
		[p[1][2] for p in poincare_section], 
		[p[1][4] for p in poincare_section],
		marker_z=[p[2] for p in poincare_section],
		markerstrokewidth=0.0,
		markersize=1,
		mc=:blues,
		cbar_title="\$ln(t)\$ for Lyapunov spectrum convergance",
		label=false,
		xlabel="\$y\$", ylabel="\$v_y\$",
		title="Poincaré section of Hénon–Heiles system",
		backgroundinside=:linen
	)
end
savefig("plots/exercise3_15.png")
