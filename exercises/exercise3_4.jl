using ForwardDiff, LinearAlgebra, Plots

abstract type DynamicalSystem end

mutable struct DiscreteDynamicalSystem <: DynamicalSystem
	rule::Function
	u::Vector{Float64}
end

mutable struct ContinuousDynamicalSystem <: DynamicalSystem
	rule::Function
	u::Vector{Float64}
end

function lyapunovspectrum(ds::DiscreteDynamicalSystem, N::Int64, Δt::Int64, k::Int64)
	# N = 10_000, Δt = 5, Δt << N & Δt | N
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

#=
function sm_rule(u)
	θ, v = u
	k = 0.6
	return [
		θ + v + k * sin(θ),
		    v + k * sin(θ)
	]
end
sm = DiscreteDynamicalSystem(sm_rule, [0.1, 0.11])
lyapunovspectrum(sm, 10_000, 5, 2)
=#

λ1s = []
λ2s = []
for k=0.6:0.01:1.2
	function sm_rule_k(u)
		θ, v = u
		return [
			θ + v + k * sin(θ),
			    v + k * sin(θ)
		]
	end

	ds = DiscreteDynamicalSystem(sm_rule_k, [0.1, 0.11])
	λ1, λ2 = lyapunovspectrum(ds, 10_000, 5, 2)
	push!(λ1s, λ1)
	push!(λ2s, λ2)
end

plot(
	0.6:0.01:1.2, λ1s,
	xlabel="k", ylabel="λ", label="λ1",
	lc=:red,
	title="Lyapunov Spectrum of Standard Map"
)
plot!(
	0.6:0.01:1.2, λ2s,
	xlabel="k", ylabel="λ", label="λ2",
	lc=:blue
)
plot!(
	0.6:0.01:1.2, λ1s .+ λ2s,
	xlabel="k", ylabel="λ", label="λ1 + λ2",
	lc=:orange
)
# λ1 + λ2 = 0 since the Standard Map is a conservative system
savefig("plots/exercise3_4.png")
