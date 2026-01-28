using DifferentialEquations
using Plots
using DynamicalSystems

function rössler!(du, u, p, t)
	a, b, c = p
	x, y, z = u

	du[1] = -y - z
	du[2] = x + a * y
	du[3] = b + z * (x - c)
end

prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0,1_000.0), [0.2, 0.2, 5.7])
solv = solve(prob; dtmax=0.01)

p1 = plot(solv, idxs=(1,2,3), lw=0.5, lα=0.1, label=false, title="Rössler system")

δts = range(0.1, 1.0, 10)
colors = [:silver, :red, :orange, :yellow, :green, :cyan, :blue, :purple, :pink, :brown]
p2  = plot()

for (δt, c) in zip(δts, colors)
	local solv = solve(prob; dtmax=δt)
	SMIs = []
	for τ in 1:20
		smi = selfmutualinfo(solv.u, [τ, τ, τ])[1]
		push!(SMIs, smi)
	end
	plot!(1:20, SMIs ./ maximum(SMIs), label="\$ δt = $δt\$", lc=c, lw=2, xlabel="\$τ\$", ylabel="SMI")
end
plot(p1, p2)
savefig("plots/exercise6_4.png")
