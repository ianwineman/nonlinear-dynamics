using CSV, DataFrames
using Plots
using DynamicalSystems

"Permutation of relative amplitude"
function permra(x)
	return [findfirst(==(y), sort(x)) for y in x]
end

"Permutation entropy"
function perment(X, d)
	pras = [permra(X[i:(i+d-1)]) for i in 1:(length(X)-d+1)]
	upras = unique(pras)
	Pis = [count(==(up), upras) |> float for up in upras]
	Pis ./= length(Pis)
	return -sum(Pis .* log2.(Pis))
end


df = CSV.read("book/exercise_data/7.csv", DataFrame; header=false)
X = df[:, 1]
my_PEs    = [perment(X, d) for d in 2:10]
their_PEs = [entropy_permutation(X; m=d) for d in 2:10]

plot(
	2:10, 
	my_PEs,
	xlabel="\$ d\$",
	ylabel="Permutation entropy",
	label="Mine",
	xlim=[2, 10],
	xtick=2:10,
	ylim=[0, ceil(maximum(PEs))],
	lc=RGB(159/255, 178/255, 97/255),
	lw=1.5
)
plot!(
	2:10,
	their_PEs,
	label="DynamicalSystems.jl",
	lc=RGB(96/255, 77/255, 158/255),
	lw=1.5
)
savefig("plots/exercise6_13.png")
