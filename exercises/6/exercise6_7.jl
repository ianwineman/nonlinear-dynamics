using LinearAlgebra
using DynamicalSystems
using CSV, DataFrames
using Plots

function nearest_neighbor(x,X; w=5)
	i = findfirst(==(x), X)
	norms = [norm(x - y) for y in X]
	norms[maximum([1, i-w]):minimum([length(X), i+w])] .= Inf
	return argmin(norms)
end

df = CSV.read("book/exercise_data/2.csv", DataFrame; header=false)
R = reduce(hcat, eachcol(df)) |> 
	eachrow |> 
	collect |> 
	x -> map(y -> Vector(y), x)

function lyapunov_from_data(R, k, δt)
	M = length(R) - k
	ξs = zeros(k)

	Ζ = [(i,nearest_neighbor(R[i],R[1:M])) for i in 1:M] |>
		x -> map(y -> [y .+ (i,i) for i=0:k-1], x) |> 
		stack |>
		x -> map(y -> norm(R[y[1]] - R[y[2]]), x)

	for n = 1:k
		ξn = sum([log(ζ) for ζ in Ζ[n, :]]) / M
		ξs[n] = ξn
	end
	return ξs ./ δt
end

ks = 1:150
lfd = lyapunov_from_data(R, last(ks), 0.05)
λ = (lfd[120]-lfd[20])/(ks[120]-ks[20])

ds_lfd = DynamicalSystems.lyapunov_from_data(StateSpaceSet(R), ks; w=5)
ds_λ = ChaosTools.linreg(collect(ks) .* 0.05, ds_lfd)[2]

plot(
	ks, 
	lfd, 
	label="\$ξ_n/δt\$ versus \$n\$",
	xlabel="\$n\$",
	ylabel="\$ξ_n/δt\$",
	lw=2,
	lc=:blue,
	xtick=0:10:150,
	xlim=[0,150],
	title="Largest Lyapunov exponent from data"
)
plot!(
	1:150, 
	collect(1:150) .*    λ .+ lfd[1], 
	label="\$ λ ≈ $(round(λ; digits=3))\$ (mine)",
	lc=:red
)
plot!(
	1:150, 
	collect(1:150) .* ds_λ .+ lfd[1], 
	label="\$ λ ≈ $(round(ds_λ; digits=3))\$ (DynamicalSystems.jl)",
	lc=:orange
)
savefig("plots/exercise6_7.png")
