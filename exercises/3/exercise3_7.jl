using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots 
using PredefinedDynamicalSystems, DynamicalSystems

abstract type MyDynamicalSystem end

mutable struct MyContinuousDynamicalSystem <: MyDynamicalSystem
	rule::Function     # f(u) (works for ForwardDiff)
	rule_de::Function  # f(du,u,p,t) (works for DifferentialEquations) 
	u::Vector{Float64} #
	p::Vector{Float64} # for DifferentialEquations
end

function Mylyapunovspectrum(ds::MyContinuousDynamicalSystem, N::Float64, Δt::Float64, k::Int64)
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

Myλs = []
λs = []
for ρ=20:100
	function lorenz63_rule(u)
		x, y, z = u
		σ, ρ, β = [10.0, ρ, 8/3]
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

	ds = MyContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, [20.0, 20.0, 20.0], [10.0, ρ, 8/3])
	push!(Myλs, Mylyapunovspectrum(ds, 1_000.0, 1.0, 3))

	lo = PredefinedDynamicalSystems.lorenz([20.0, 20.0, 20.0]; σ=10.0, ρ=ρ, β=8/3)
	push!(λs, lyapunovspectrum(lo, 1_000))
end

plot(
	collect(20:100), 
	hcat([λ[1] for λ in Myλs], [λ[2] for λ in Myλs], [λ[3] for λ in Myλs], [λ[1] for λ in λs], [λ[2] for λ in λs], [λ[3] for λ in λs]),
	xlabel="ρ", ylabel="λ", label=["λ1" "λ2" "λ3" "λ1" "λ2" "λ3"],
	lc=[:red :orange :blue :red :orange :blue],
	ls=[:solid :solid :solid :dot :dot :dot],
	legendtitle="Solid: mine\nDotted: DynamicalSystems.jl",
	legendtitlefontsize=4, legendfontsize=4,
	legendposition=:outerbottomright,
	title="Lyapunov Spectrum of Lorenz63 system"
)
savefig("plots/exercise3_7.png")
