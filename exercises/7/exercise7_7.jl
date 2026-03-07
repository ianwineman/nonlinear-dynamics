using CSV, DataFrames
using Statistics
using Associations

for file in [14, 15]
	df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
	X, Y = df[:, 1], df[:, 2]

	xtoy_tes = []
	ytox_tes = []
	for delay in 1:10
		precise = true
		discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise)))
		est_disc_to = EntropyDecomposition(TEShannon(; embedding=EmbeddingTE(; dS=delay, dT=delay)), PlugIn(Shannon()), discretization);
		xtoy_te = association(est_disc_to, X, Y)
		ytox_te = association(est_disc_to, Y, X)
		push!(xtoy_tes, xtoy_te)
		push!(ytox_tes, ytox_te)
	end
	println("$file.csv, x->y $(mean(xtoy_tes)), y->x $(mean(ytox_tes))")
end

# 14.csv, x->y 0.02571056416239622, y->x 0.02670107285692843
# 15.csv, x->y 0.002426487682476172, y->x 0.0019092320588688233
