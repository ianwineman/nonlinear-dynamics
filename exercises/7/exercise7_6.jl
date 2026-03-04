using Associations
using DifferentialEquations
using Plots, StatsPlots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end

prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,20.0), [10.0, 40.0, 8/3])
solv = solve(prob)

ϵs = [0.1, 0.5, 1.0]
transferentropies = []
for (i, j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
	tes = []
	for ϵ in ϵs
		println("$i -> $j")
		precise = true
		discretization = CodifyVariables(TransferOperator(RectangularBinning(ϵ, precise)))
		est_disc_to = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
		te = association(est_disc_to, getindex.(solv.u, i), getindex.(solv.u, j))
		push!(tes, te)
	end
	push!(transferentropies, tes)
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
p2 = groupedbar(
	["\$ $(vars[i]) → $(vars[j])\$" for (i,j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]],
	Float64.(stack(transferentropies; dims=1)),
	label="\$r=" .* hcat(string.(ϵs)...) .* "\$",
	xlabel="Transfer direction",
	ylabel="Transfer entropy"
)

plot(p1, p2)
savefig("plots/exercise7_6.png")
