using CSV, DataFrames
using LinearAlgebra
using Plots

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

files = [2, 3, 4, 6, 9]
Δs = zeros(length(files))

for (i,f) in enumerate(files)
	println("$i/$(length(files))")
	df = CSV.read("book/exercise_data/$f.csv", DataFrame; header=false)
	X = [collect(df[i, :]) for i in 1:size(df)[1]]
	Δs[i] = correlation_dimension(X, 3)
end

scatter(
	files, Δs, 
	xticks=[2,3,4,6,9], 
	xlabel="Dataset number", 
	yticks=round.(Δs; digits=2),
	ylabel="Fractal dimension", 
	label=false, 
	msw=0, ms=5, mc=:hotpink
)
savefig("plots/exercise5_11.png")
