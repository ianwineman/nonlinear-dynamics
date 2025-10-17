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


let u0 = [20.0, 20.0, 20.0], p  = [10.0, 28.0, 8/3]
	function lorenz63_rule(u)
		x, y, z = u
		σ, ρ, β = [10.0, 28.0, 8/3]
		return [
			σ * (y - x),
			-x * z + ρ * x - y,
			x * y - β * z
		]
	end

	function lorenz63_rule_de!(du, u, p, t)
		x, y, z = u
		σ, ρ, β = p
		du[1] = σ * (y - x)
		du[2] = -x * z + ρ * x - y
		du[3] = x * y - β * z
	end

	lo = ContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, u0, p)
	lsN = lyapunovspectrum_convergance(lo, 1_000.0, 1.0, 3)

	plot(
		1:length(lsN), hcat([ls[1] for ls in lsN], [ls[2] for ls in lsN], [ls[3] for ls in lsN]),
		xlabel="Iterations, \$N\$", ylabel="\$λ\$",
		label=hcat(
			"\$λ_1 = $(round(lsN[end][1], digits=4))\$",
			"\$λ_2 = $(round(lsN[end][2], digits=4))\$",
			"\$λ_3 = $(round(lsN[end][3], digits=4))\$",
		),
		linecolor=[:red :orange :blue],
		title="Converagance of Lyapunov Spectrum"
	)
	scatter!(
		[length(lsN), length(lsN), length(lsN)], [lsN[end][1], lsN[end][2], lsN[end][3]], 
		label=false,
		markercolor=[:red, :orange, :blue]
	)
	savefig("plots/exercise3_14_1.png")
end

let Δts = [0.1, 0.5, 1.0, 5.0, 10.0]
	function lorenz63_rule(u)
		x, y, z = u
		σ, ρ, β = [10.0, 28.0, 8/3]
		return [
			σ * (y - x),
			-x * z + ρ * x - y,
			x * y - β * z
		]
	end

	function lorenz63_rule_de!(du, u, p, t)
		x, y, z = u
		σ, ρ, β = p
		du[1] = σ * (y - x)
		du[2] = -x * z + ρ * x - y
		du[3] = x * y - β * z
	end

	u0 = [20.0, 20.0, 20.0]
	p  = [10.0, 28.0, 8/3]
	lo = ContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, u0, p)

	lsNs = []
	for Δt in Δts
		lsN = lyapunovspectrum_convergance(lo, 1_000.0, Δt, 3)
		push!(lsNs, lsN)
	end
	plot(
		1:length(lsNs[1]),
		hcat(
			[[ls[1] for ls in lsN] for lsN in lsNs]..., 
			[[ls[2] for ls in lsN] for lsN in lsNs]..., 
			[[ls[3] for ls in lsN] for lsN in lsNs]...
		),
		xlabel="Iterations, \$N\$", ylabel="\$λ\$", 
		label=false,
		#linecolor=[reds(n), oragnes(n), blues(n)],
		title="Converagance of Lyapunov Spectrum",
		legendposition=:outertopright,
		legendtitle="\$Δt∈$Δts\$"
	)
	scatter!(
		[length(lsNs[1])], [lsNs[3][end][1]],
		label="\$λ_1 = $(round(lsNs[3][end][1], digits=4))\$",
		markercolor=:red
	)
	scatter!(
		[length(lsNs[1])], [lsNs[3][end][2]],
		label="\$λ_2 = $(round(lsNs[3][end][2], digits=4))\$",
		markercolor=:orange
	)
	scatter!(
		[length(lsNs[1])], [lsNs[3][end][3]],
		label="\$λ_3 = $(round(lsNs[3][end][3], digits=4))\$",
		markercolor=:blue
	)
	savefig("plots/exercise3_14_2.png")
end
