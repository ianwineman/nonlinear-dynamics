using LinearAlgebra
using DynamicalSystems
using CSV, DataFrames
using Plots, Measures

function nearest_neighbor(x,X; w=5)
	i = findfirst(==(x), X)
	norms = [norm(x - y) for y in X]
	norms[maximum([1, i-w]):minimum([length(X), i+w])] .= Inf
	return argmin(norms)
end

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

### File: 3.csv; δt = 0.01 ###
df = CSV.read("book/exercise_data/3.csv", DataFrame; header=false)
R = reduce(hcat, eachcol(df)) |> 
	eachrow |> 
	collect |> 
	x -> map(y -> Vector(y), x)

ks = 1:150
lfd = lyapunov_from_data(R, last(ks), 0.01)
λ = (lfd[120]-lfd[40])/(ks[120]-ks[40])

ds_lfd = DynamicalSystems.lyapunov_from_data(StateSpaceSet(R), ks; w=5)
ds_λ = ChaosTools.linreg(collect(ks) .* 0.01, ds_lfd)[2]

p1 = plot(
	ks, 
	lfd, 
	label="\$ξ_n/δt\$ versus \$n\$ (3.csv)",
	xlabel="\$n\$",
	ylabel="\$ξ_n/δt\$",
	lw=2,
	lc=:blue,
	xtick=0:20:150,
	xlim=[1,150],
	left_margin=10mm,
	bottom_margin=5mm
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

### File: 6.csv; δt = 1 ###
df = CSV.read("book/exercise_data/6.csv", DataFrame; header=false)
R = reduce(hcat, eachcol(df)) |> 
	eachrow |> 
	collect |> 
	x -> map(y -> Vector(y), x)

ks = 1:15
lfd = lyapunov_from_data(R, last(ks), 1.0)
λ = (lfd[10]-lfd[1])/(ks[10]-ks[1])

ds_lfd = DynamicalSystems.lyapunov_from_data(StateSpaceSet(R), ks; w=5)
ds_λ = ChaosTools.linreg(collect(ks)[1:10] .* 1.0, ds_lfd[1:10])[2]

p2 = plot(
	ks, 
	lfd, 
	label="\$ξ_n/δt\$ versus \$n\$ (6.csv)",
	xlabel="\$n\$",
	ylabel="\$ξ_n/δt\$",
	lw=2,
	lc=:blue,
	xtick=0:15,
	xlim=[1,15],
	left_margin=10mm,
	bottom_margin=5mm
)
plot!(
	1:15, 
	collect(1:15) .*    λ .+ lfd[1], 
	label="\$ λ ≈ $(round(λ; digits=3))\$ (mine)",
	lc=:red
)
plot!(
	1:15, 
	collect(1:15) .* ds_λ .+ lfd[1], 
	label="\$ λ ≈ $(round(ds_λ; digits=3))\$ (DynamicalSystems.jl)",
	lc=:orange
)

### File: 10.csv; δt = 0.01 ###
df = CSV.read("book/exercise_data/10.csv", DataFrame; header=false)
R = reduce(hcat, eachcol(df)) |> 
	eachrow |> 
	collect |> 
	x -> map(y -> Vector(y), x)

ks = 1:150
lfd = lyapunov_from_data(R, last(ks), 0.01)
λ = (lfd[120]-lfd[30])/(ks[120]-ks[30])

ds_lfd = DynamicalSystems.lyapunov_from_data(StateSpaceSet(R), ks; w=5)
ds_λ = ChaosTools.linreg(collect(ks)[30:120] .* 0.01, ds_lfd[30:120])[2]

p3 = plot(
	ks, 
	lfd, 
	label="\$ξ_n/δt\$ versus \$n\$ (10.csv)",
	xlabel="\$n\$",
	ylabel="\$ξ_n/δt\$",
	lw=2,
	lc=:blue,
	xtick=0:20:150,
	xlim=[1,150],
	left_margin=10mm,
	bottom_margin=5mm
)
plot!(
	30:120, 
	collect(0:90) .*    λ .+ lfd[30], 
	label="\$ λ ≈ $(round(λ; digits=3))\$ (mine)",
	lc=:red
)
plot!(
	30:120, 
	collect(0:90) .* ds_λ .+ lfd[30], 
	label="\$ λ ≈ $(round(ds_λ; digits=3))\$ (DynamicalSystems.jl)",
	lc=:orange
)

plot(
	p1, p2, p3, 
	layout=(1,3), 
	size=(1200,400),
	plot_title="Largest Lyapunov exponent from data"
)
savefig("plots/exercise6_8.png")
