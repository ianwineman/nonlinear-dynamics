using PredefinedDynamicalSystems
using DifferentialEquations, LinearAlgebra, ForwardDiff
using Plots

abstract type DynamicalSystem end

mutable struct ContinuousDynamicalSystem <: DynamicalSystem
	rule::Function     # f(u) (works for ForwardDiff)
	rule_de::Function  # f(du,u,p,t) (works for DifferentialEquations) 
	u::Vector{Float64} #
	p::Vector{Float64} # for DifferentialEquations
end

function correlation_sum(X, ϵ, w)
	N = length(X)
	C = 0

	for i in 1:(N - w - 1)
		for j in (1 + w + i):N
			C += ifelse(norm(X[i] .- X[j]) < ϵ, 1, 0)
		end
	end

	return (2*C)/((N - w)*(N - w - 1))
end

function correlation_dimension(X, w)
	ϵs = [exp(j) for j in -10.0:0.1:5.0]
	Cs = [correlation_sum(X, ϵ, w) for ϵ in ϵs]

	i = findfirst(x -> !isinf(x), log.(Cs)) + round(Int, length(Cs)/10)
	j = findfirst(x -> x == 0.0, log.(Cs))  - round(Int, length(Cs)/10)
	Δ = (log.(Cs)[j] - log.(Cs)[i]) / (log.(ϵs)[j] - log.(ϵs)[i])

	return Δ
end

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

function lyapunovdimension(λs)
	sums = [sum(λs[1:j]) for j in 1:length(λs)]
	k = findfirst(s->s<0, sums)
	if k == 1 || k == nothing
		return 0.0
	end
	k -= 1
	return k + sums[k]/abs(λs[k+1])
end

function henonheiles!(du, u, p, t)
	x, y, vx, vy = u
	du[1] = vx
	du[2] = vy
	du[3] = -x - 2*x*y
	du[4] = -y - (x^2 - y^2)
end

u0s = first(PredefinedDynamicalSystems.henonheiles_ics(0.13, 5), 5)
short_trajectories = []
full_trajectories = []
for u0 in u0s
	prob = ODEProblem(henonheiles!, u0, (0.0,1_000.0))
	solv = solve(prob; dtmax=0.1)
	trajectory = first(solv.u, 10_000)
	push!(full_trajectories, trajectory)
	for i = 1:10
		push!(short_trajectories, trajectory[((i-1)*1_000 + 1):(i*1_000)])
	end
end

dimensions = Vector{Tuple{Float64, Float64}}(undef, length(short_trajectories)+length(full_trajectories))
Threads.@threads for (i, tr) in collect(enumerate(vcat(short_trajectories, full_trajectories)))
	println("start $i on thread$(Threads.threadid())")

	ΔC = correlation_dimension(tr, 3)
	hh = ContinuousDynamicalSystem(henon_heiles_rule, henon_heiles_rule_de!, tr[1], [])
	ls = lyapunovspectrum(hh, Float64(length(tr))/10, 1.0, 4)
	ΔL = lyapunovdimension(ls)

	dimensions[i] = (ΔC, ΔL)
end

colors = vcat(
	[:red for _ in 1:10], 
	[:orange for _ in 1:10], 
	[:blue for _ in 1:10], 
	[:green for _ in 1:10], 
	[:purple for _ in 1:10], 
	[:red, :orange, :blue, :green, :purple]
)
scatter(
	first.(enumerate(dimensions)), 
	getindex.(last.(enumerate(dimensions)), 1), 
	label="ΔC (color denotes u0)", 
	xtick=false, 
	ylabel="Dimension", 
	xlabel="Trajectory", 
	title="Hénon-Heiles system", 
	mc=colors
)
scatter!(
	first.(enumerate(dimensions)), 
	getindex.(last.(enumerate(dimensions)), 2), 
	label="ΔL (color denotes u0)", 
	shape=:rect, 
	mc=colors
)
vline!(
	[50.5], 
	label="length = 1,000 | length = 10,000", 
	style=:dash, 
	lc=:black
)
savefig("plots/exercise5_13.png")
