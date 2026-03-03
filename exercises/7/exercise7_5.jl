using Associations
using DifferentialEquations
using Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end

prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,20.0), [10.0, 40.0, 8/3])
solv = solve(prob)

transferentropies = []
for (i, j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
	println("$i -> $j")
	precise = true
	discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise)))
	est_disc_to = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
	te = association(est_disc_to, getindex.(solv.u, i), getindex.(solv.u, j))
	push!(transferentropies, te)
end

p1 = plot(
	solv, 
	idxs=(1,2,3),
	#label="Lorenz63",
	label="\$ ẋ = σ(y-x)\$ \n\$ ẏ = -xz+ρx-y\$ \n\$ ż = xy - βz\$",
	xlabel="\$ x\$",
	ylabel="\$ y\$",
	zlabel="\$ z\$"
)

vars = ["x", "y", "z"]
p2 = bar(
	["\$ $(vars[i]) → $(vars[j])\$" for (i,j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]],
	transferentropies,
	label=false,
	xlabel="Transfer direction",
	ylabel="Transfer entropy"
)

plot(p1, p2)
savefig("plots/exercise7_5.png")
