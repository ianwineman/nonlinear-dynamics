using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots

abstract type DynamicalSystem end

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


function lyapunovdimension(λs)
	sums = [sum(λs[1:j]) for j in 1:length(λs)]
	k = findfirst(s->s<0, sums) - 1
	if k == 0
		return 0.0
	end
	return k + sums[k]/abs(λs[k+1])
end


λs = []
Δs = []
for ρ=20:100
	println(ρ)
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

	ds = ContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, [20.0, 20.0, 20.0], [10.0, ρ, 8/3])
	ls = lyapunovspectrum(ds, 1_000.0, 1.0, 3)
	push!(λs, ls)
	push!(Δs, lyapunovdimension(ls))
end

plot(
	collect(20:100), 
	hcat([λ[1] for λ in λs], [λ[2] for λ in λs], [λ[3] for λ in λs]),
	ytick=collect(-18:4:2),
	xlabel="ρ", ylabel="λ", label=["λ1" "λ2" "λ3"],
	lc=[:red :orange :blue],
	title="Lyapunov spectrum & dimension (Lorenz63)",
	legendposition=:right
)
plot!(
	collect(20:100), Δs, 
	lc=:mediumpurple2, 
	lw=2,
	label="ΔL"
)
savefig("plots/exercise5_12.png")


function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end
prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,100.0), [10.0, 40.0, 8/3])
solv = solve(prob)
plot(solv, idxs=(1,2,3), label="ρ=40, ΔL ≈ 2")

prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,100.0), [10.0, 22.0, 8/3])
solv = solve(prob)
plot!(solv, idxs=(1,2,3), label="ρ=22, ΔL ≈ 0")
savefig("plots/exercise5_12_2.png")
