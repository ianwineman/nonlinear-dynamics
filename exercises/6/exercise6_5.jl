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

files = [1, 3, 5, 12]

Xs = []
Ys = []

for file in files
	df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
	X = df[:, 1]
	push!(Xs, X)
	Y = []
	for τ in 0:50
		push!(Y, selfmutualinfo(X, [τ, τ, τ])[1])
	end
	push!(Ys, Y ./ maximum(Y))
end

p1 = plot(
	0:50,
	Ys,
	xtick=0:2:50, 
	xrotation=45, 
	xtickfontsize=4, 
	xlim=[0,50], 
	ylim=[0,1],
	xguide="\$τ\$",
	yguide="\$ SMI\$",
	label=["\$1\$" "\$3\$" "\$5\$" "\$12\$"],
	lw=2,
	lc=[:red :orange :blue :purple]
)

τs = [55, 17, 20, 12] # from p1

Eds = []

for (i, w) in enumerate(Xs)
	push!(Eds, delay_afnn(w, 1, 1:10; w=5))
end

p2 = plot(
	1:10,
	Eds,
	xlabel="\$ d\$", 
	ylabel="\$E_d\$",
	xtick=1:10,
	label=["\$1\$" "\$3\$" "\$5\$" "\$12\$"],
	lw=2,
	lc=[:red :orange :blue :purple]
)

plot(p1, p2, size=(800,400))
savefig("plots/exercise6_5.png")
