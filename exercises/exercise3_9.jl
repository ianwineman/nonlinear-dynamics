using DifferentialEquations, ForwardDiff, LinearAlgebra, Plots

abstract type DynamicalSystem end

mutable struct DiscreteDynamicalSystem <: DynamicalSystem
	rule::Function
	u::Vector{Float64}
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

function towel_rule(u)
	x, y, z = u
	return [
		3.8*x*(1.0-x) - 0.05*(y+0.35)*(1.0-2.0*z),
		0.1*((y+0.35)*(1.0-2.0*z)-1.0)*(1.0-1.9*x),
		3.78*z*(1.0-z) + 0.2*y
	]
end
towel = DiscreteDynamicalSystem(towel_rule, [0.085, -0.121, 0.075])
lyapunovspectrum(towel, 10_000, 5, 3)
