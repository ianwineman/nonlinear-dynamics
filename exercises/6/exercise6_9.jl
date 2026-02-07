using LinearAlgebra
using CSV, DataFrames
using Plots, Measures
using ProgressMeter

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

function correlation_dim(X; ϵs=[exp(i) for i in -10.0:0.1:5.0], w=3)
	Cs = []
	for ϵ in ϵs
		push!(Cs, correlation_sum(X, ϵ, w))
	end
	i = findfirst(x -> !isinf(x), log.(Cs)) + round(Int, length(Cs)/10)
	j = findfirst(x -> x == 0.0, log.(Cs))  - round(Int, length(Cs)/10)
	xstart = ϵs[i]
	xend = ϵs[j]
	Δ = (log.(Cs)[j] - log.(Cs)[i]) / (log.(ϵs)[j] - log.(ϵs)[i])
	return Δ
end

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

files = [3, 6, 12]
colors = [:red :orange :blue]

SMIs = []
for file in files
	df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
	R = reduce(hcat, eachcol(df)) |> 
		eachrow |> 
		collect |> 
		x -> map(y -> Vector(y), x)

	SMI = selfmutualinfo(R, 1:30)
	push!(SMIs, SMI)
end

p1 = plot(
	1:30, 
	SMIs, 
	xtick=1:30,
	xtickfontsize=7,
	xlim=[1,30],
	xlabel="\$τ\$",
	ylabel="\$SMI\$",
	labels=string.(files'),
	legendtitle="File",
	lw=2,
	lc=colors,
	left_margin=5mm,
	bottom_margin=5mm
)

τs = [16, 10, 12] # from p1

Δss = []
for (i,file) in enumerate(files)
	df = CSV.read("book/exercise_data/$file.csv", DataFrame; header=false)
	R = reduce(hcat, eachcol(df)) |> 
		eachrow |> 
		collect |> 
		x -> map(y -> Vector(y), x)

	Δs = zeros(10)
	@showprogress desc="[$file]" Threads.@threads for d in 1:10
		emb = delay_embed(first.(R), τs[i], d)
		Δ = correlation_dim(emb; ϵs=[exp(i) for i in -10.0:0.3:5.0])
		Δs[d] = Δ
	end
	push!(Δss, Δs)
end

p2 = plot(
	1:10,
	Δss,
	xtick=1:10,
	xlim=[1,10],
	xlabel="Embedding dimension \$ d\$",
	ylabel="Fractal dimension \$ Δ\$",
	labels=string.(files') .* reshape([", \$ τ=$t\$" for t in τs], (1,3)),
	legendtitle="File",
	lw=2,
	lc=colors,
	left_margin=5mm,
	bottom_margin=5mm
)

plot(p1, p2, size=(1000,400))
savefig("plots/exercise6_9.png")
