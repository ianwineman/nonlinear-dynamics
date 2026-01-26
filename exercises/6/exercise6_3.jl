using Plots
using Statistics: mean
using LinearAlgebra: norm
using DynamicalSystems

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

function nn(vi, vis; w=5)
	i = findfirst(==(vi), vis)
	window = (i-w):(i+w)
	norms = [norm(vis[j] - vi) for j in 1:length(vis)]
	n1 = norms[1:first(window)-1]
	n2 = norms[last(window)+1:end]
	if minimum(n1; init=typemax(Int64)) >= minimum(n2; init=typemax(Int64))
		return argmin(n2)
	else
		return argmin(n1)
	end
end

function my_delay_afnn(w, τ::Int64, ds::UnitRange{Int64})
	ads = []
	for d in first(ds):(last(ds)+1)
		vis_d    = delay_embed(w, τ, d)
		vinns_d  = Dict([i => nn(vi, vis_d) for (i, vi) in enumerate(vis_d)])

		vis_d1   = delay_embed(w, τ, d+1)
		vinns_d1 = [vis_d1[vinns_d[i]] for i in 1:length(vinns_d) if vinns_d[i] <= length(vis_d1)]

		ad = [norm(vis_d1[i] - vinns_d1[i])/norm(vis_d[i] - vis_d[vinns_d[i]]) for i in minimum(length.([vis_d, vinns_d, vis_d1, vinns_d1]))]
		push!(ads, mean(ad))
	end
	return [ads[i+1]/ads[i] for i in 1:(length(ads)-1)]
end

function towel_map(u)
	x, y, z = u
	return [
		3.8 * x * (1 - x) - 0.05 * (y + 0.35) * (1 - 2 * z),
		0.1 * ((y + 0.35) * (1 - 2 * z) - 1) * (1 - 1.9 * x),
		3.78 * z * (1 - z) + 0.2 * y
	]
end

U = [[0.085, -0.121, 0.075]]
for _ in 1:1_000
	push!(U, towel_map(U[end]))
end

w = first.(U)

Eds = my_delay_afnn(w, 1, 1:10)
p1 = plot(1:10, Eds, label="my_delay_afnn", xlabel="\$τ\$", ylabel="\$E_d\$", xtick=1:10)
plot!(1:10, delay_afnn(w, 1, 1:10; w=5), label="DynamicalSystems.jl")

emb = delay_embed(w, 1, 3)
p2 = plot(
	getindex.(emb, 1),
	getindex.(emb, 2),
	getindex.(emb, 3),
	lα=0.1,
	camera=(135,30),
	m=:circ,
	axis=false,
	grid=false,
	label=false,
	msw=0.2,
	title="Towel map, \$x\$ delay embedding"
)
savefig("plots/exercise6_3.png")
