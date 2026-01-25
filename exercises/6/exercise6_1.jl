using DifferentialEquations
using Plots, Plots.PlotMeasures, LaTeXStrings

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

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end

prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,20.0), [10.0, 40.0, 8/3])
solv = solve(prob; dtmax=0.01)

p1 = plot(
	solv, 
	idxs=(1,2,3), 
	label=false, 
	lw=0.5,
	lc=:black,
	showaxis=false,
	title="Lorenz63\nDelay Embeddings",
	titlelocation=:left
)

## (i)
τ = [0, 17, 34]
j = [1, 1, 1]

emb = delay_embed(solv.u, τ, j)

p2 = plot(
	getindex.(emb, 1), 
	getindex.(emb, 2), 
	getindex.(emb, 3), 
	title=L"(i)\, \textbf{j} = (1,1,1)", 
	label=false, 
	titlefontsize=8,
	lw=0.5,
	lc=:red,
	showaxis=false
)

## (ii)
τ = [0, 15, 30]
j = [3, 3, 3]

emb = delay_embed(solv.u, τ, j)

p3 = plot(
	getindex.(emb, 1), 
	getindex.(emb, 2), 
	getindex.(emb, 3), 
	title=L"(ii)\, \textbf{j} = (3,3,3)", 
	label=false, 
	titlefontsize=8,
	lw=0.5,
	lc=:orange,
	showaxis=false
)

## (iii)
τ = [0, 17, 15]
j = [1, 1, 3]

emb = delay_embed(solv.u, τ, j)

p4 = plot(
	getindex.(emb, 1), 
	getindex.(emb, 2), 
	getindex.(emb, 3), 
	title=L"(iii)\, \textbf{j} = (1,1,3)", 
	label=false, 
	titlefontsize=8,
	lw=0.5,
	lc=:blue,
	showaxis=false
)

## (iv)
τ = [0, 15, 30]
j = [2, 3, 3]

emb = delay_embed(solv.u, τ, j)

p5 = plot(
	getindex.(emb, 1), 
	getindex.(emb, 2), 
	getindex.(emb, 3), 
	title=L"(iv)\, \textbf{j} = (2,3,3)", 
	label=false, 
	titlefontsize=8,
	lw=0.5,
	lc=:purple,
	showaxis=false
)

##
plot(
	p1, p2, p3, p4, p5, 
	layout=@layout(
		[a{0.3w} [grid(2,2)]]
	)
)
savefig("plots/exercise6_1.png")
