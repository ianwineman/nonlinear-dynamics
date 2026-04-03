using DifferentialEquations
using Plots, StatsPlots
using LinearAlgebra
using Associations

function lorenz96!(du, u, p, t)
	D, F = Int(p[1]), p[2]

    for i in 1:D
    	du[i] = ((u[mod1(i+1, D)] - u[mod1(i-2, D)]) * u[mod1(i-1, D)]) - u[i] + F
    end
end

u0 = [0.12013623147677799, 0.274707399757551, 0.4738961112067477, 0.14933883808016557, 0.4274098963989116] # rand(5)
prob = ODEProblem(lorenz96!, u0, (0.0, 100.0), [5.0, 8.0])
solv = solve(prob)

rs = 2:5
transferentropies = []
for (i, j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]
	tes = []
	for r in rs
		println("$i -> $j")
		precise = true
		discretization = CodifyVariables(TransferOperator(RectangularBinning(r, precise)))
		est_disc_to = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
		te = association(est_disc_to, getindex.(solv.u, i), getindex.(solv.u, j))
		push!(tes, te)
	end
	push!(transferentropies, tes)
end

vars = ["x_1", "x_2", "x_3"]
groupedbar(
	["\$ $(vars[i]) → $(vars[j])\$" for (i,j) in [(1, 2), (1, 3), (2, 1), (2, 3), (3, 1), (3, 2)]],
	Float64.(stack(transferentropies; dims=1)),
	label="\$r=" .* hcat(string.(rs)...) .* "\$",
	xlabel="Transfer direction",
	ylabel="Transfer entropy",
	title="Lorenz96 \$ (D=5, F=8)\$"
)
savefig("plots/exercise7_14.png")
