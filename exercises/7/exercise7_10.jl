using CSV, DataFrames
using Statistics, LinearAlgebra
using DynamicalSystems
using Plots

function delayembed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

function nearestneighbors(x, X; k=1, w=5)
	i = findfirst(==(x), X)
	if i == nothing
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
	else
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		window = maximum([1, i-w]):minimum([length(X), i+w])
		norms[window] .= [(typemax(Int), Inf) for _ in window]
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
	end
end

function convergentcrossmapping(x, y, τ, d)
	mx = delayembed(x, τ, d)
	ψ = zeros(length(y))
	for (i, x) in enumerate(mx)
		nn = nearestneighbors(x, mx; k=5)
		u = [exp(-norm(x - mx[j]) / norm(x - mx[nn[end]])) for j in nn]
		w = u ./ sum(u)
		ψ[i] = sum(w .* y[nn])
	end
	return cor(y, ψ)
end

df = CSV.read("book/exercise_data/16.csv", DataFrame; header=false, delim=" ", ignorerepeated=true)
X, Y, Z = df[:, 2], df[:, 3], df[:, 5]

τ = 1 #plot(1:30, selfmutualinfo(first.(solv.u), 1:30), xtick=1:30)
d = 4 #plot(1:10, delay_afnn(X, τ, 1:10; w=5))

plot(X, Y, Z)

ccms = []
for Ns in 500:1_000:length(X)
	ccm = convergentcrossmapping(X[1:Ns], Z[1:Ns], τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
println()
plot(
	first.(ccms), 
	last.(ccms), 
	label="\$ x → z\$",
	xlabel="\$ N_s\$",
	ylabel="CCM measure",
	legendposition=:right
)

ccms = []
for Ns in 500:1_000:length(X)
	ccm = convergentcrossmapping(Y[1:Ns], Z[1:Ns], τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
println()
plot!(
	first.(ccms), 
	last.(ccms), 
	label="\$ y → z\$",
	xlabel="\$ N_s\$",
	ylabel="CCM measure",
	legendposition=:right
)

ccms = []
for Ns in 500:1_000:length(X)
	ccm = convergentcrossmapping(Z[1:Ns], X[1:Ns], τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
println()
plot!(
	first.(ccms), 
	last.(ccms), 
	label="\$ z → x\$",
	xlabel="\$ N_s\$",
	ylabel="CCM measure",
	legendposition=:right
)

ccms = []
for Ns in 500:1_000:length(X)
	ccm = convergentcrossmapping(Z[1:Ns], Y[1:Ns], τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(
	first.(ccms), 
	last.(ccms), 
	label="\$ z → y\$",
	xlabel="\$ N_s\$",
	ylabel="CCM measure",
	legendposition=:right
)

savefig("plots/exercise7_10.png")
