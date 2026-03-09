using CSV, DataFrames
using Associations
using Plots, Measures

results = []
for dS in 1:3
	for dT in 1:3
		for bin in 2:4
			for τS in -3:-1
				for τT in -3:-1
					for file in [14, 15]
						df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
						X, Y = df[:, 1], df[:, 2]
						precise = true
						discretization = CodifyVariables(TransferOperator(RectangularBinning(bin, precise)))
						est_disc_to = EntropyDecomposition(TEShannon(; embedding=EmbeddingTE(; dS=dS, dT=dT, τS=τS, τT=τT)), PlugIn(Shannon()), discretization);
						xtoy_te = association(est_disc_to, X, Y)
						ytox_te = association(est_disc_to, Y, X)
						push!(results, (file, xtoy_te, ytox_te))
					end
				end
			end
		end
	end
end

p1 = scatter(
	[r[2] for r in results if r[1] == 14], 
	label="\$ x→y\$", 
	ms=1.5, 
	msw=0,
	mc=:blue,
	legendtitle="Dataset \$14\$",
	xlabel="Parameter combination \$ (dS, dT, τS, τT, bin)\$",
	ylabel="Transfer entropy"
)
scatter!(
	[r[3] for r in results if r[1] == 14], 
	label="\$ y→x\$", 
	ms=1.5, 
	msw=0, 
	mc=:orange
)

p2 = scatter(
	[r[2] for r in results if r[1] == 15], 
	label="\$ x→y\$", 
	ms=1.5, 
	msw=0,
	mc=:green,
	legendtitle="Dataset \$15\$",
	xlabel="Parameter combination \$ (dS, dT, τS, τT, bin)\$",
	ylabel="Transfer entropy"
)
scatter!(
	[r[3] for r in results if r[1] == 15], 
	label="\$ y→x\$", 
	ms=1.5, 
	msw=0,
	mc=:purple
)

plot(p1, p2, size=(800,400), margin=5mm)
savefig("plots/exercise7_8.png")
