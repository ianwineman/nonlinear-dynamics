using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots

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

function lyapunovspectrum(ds::DiscreteDynamicalSystem, N::Int64, Δt::Int64, k::Int64)
	# N = 10_000, Δt = 5, Δt << N & Δt | N, T = N
	D = length(ds.u)
	if k > D error("k > D not allowed.") end
	Xs = [ds.u]
	for _=1:N
		push!(Xs, ds.rule(Xs[end]))
	end
	JfXs = [ForwardDiff.jacobian(ds.rule, x) for x in Xs[2:end]]
	Y0 = Matrix{Float64}(I, D, k)
	Ys = [Y0]
	Riis = []
	for n=1:N
		Yn = JfXs[n] * Ys[n]
		Q, R = qr(Yn)
		Q = Matrix(Q)
		if n % Δt == 0
			push!(Riis, log.(abs.(diag(R))))
			Yn = Q
		end
		push!(Ys, Yn)
	end
	λs = [(1/N)*sum([R[n] for R in Riis]) for n in 1:k]
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
lo = ContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, [20.0, 20.0, 20.0], [10.0, 28.0, 8/3])
ls = lyapunovspectrum(lo, 1_000.0, 1.0, 3)

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

	ds = ContinuousDynamicalSystem(lorenz63_rule, lorenz63_rule_de!, [20.0, 20.0, 20.0], [10.0, ρ, 8/3])
	push!(λs, lyapunovspectrum(ds, 1_000.0, 1.0, 3))
end

plot(
	collect(20:100), 
	[sum(λ) for λ in λs],
	xlabel="ρ", ylabel="\$\\sum_i λ_i\$",
	title="Lyapunov Spectrum of Lorenz63 system"
)
# sum({λ_i}) < 0 since Lorenz63 system is dissipative
# should have λ_i = 0 for some i...
savefig("plots/exercise3_6.png")
