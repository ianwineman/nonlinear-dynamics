using DynamicalSystems
using CSV, DataFrames
using Plots

function delay_embed end

function delay_embed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

function delay_embed(w::Vector{Vector{Float64}}, τ::Vector{Int}, j::Vector{Int}) 
	L = length(w) - maximum(τ)
	w = reduce(vcat, w')
	return [[w[i+τ[k], j[k]] for k in 1:length(τ)] for i in 1:L]
end

files = [1, 2, 6]
Xs = []
Ys = []

for file in files
	df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
	X = df[:, 1]
	push!(Xs, X)
	Y = []
	for τ in 0:100
		push!(Y, selfmutualinfo(X, [τ, τ, τ])[1])
	end
	push!(Ys, Y ./ maximum(Y))
end

p1 = plot(
	0:100,
	Ys,
	xtick=0:2:100, 
	xrotation=45, 
	xtickfontsize=2, 
	xlim=[0,100], 
	ylim=[0,1],
	xguide="\$τ\$",
	yguide="SMI",
	label=["\$1\$" "\$2\$" "\$6\$"],
	lw=2,
	lc=[:red :orange :blue]
)

τs = [55, 28, 10] # from p1
embs = [delay_embed(Xs[i], τs[i], 3) for i in 1:length(Xs)]

plots = [p1]

for i in 1:length(embs)
	p = plot(
		getindex.(embs[i], 1),
		getindex.(embs[i], 2),
		getindex.(embs[i], 3),
		label=["\$1\$, \$τ=$(τs[i])\$", "\$2\$, \$τ=$(τs[i])\$", "\$6\$, \$τ=$(τs[i])\$"][i],
		lw=0.5,
		lc=[:red :orange :blue][i],
		tick=false,
		lα=0.5
	)
	push!(plots, p)
end
plot(plots..., layout=grid(2,2))
savefig("plots/exercise6_2.png")
