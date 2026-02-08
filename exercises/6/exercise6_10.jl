using Plots
using DynamicalSystems
using Statistics
using LinearAlgebra

function henon(u; p=[1.4, 0.3])
	x, y = u
	a, b = p
	return [1 - a*x^2 + y, b*x]
end

function delay_embed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

function nearest_neighbors(x, X; k=1, w=5)
	i = findfirst(==(x), X)
	if i == nothing
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		return first.(sort(norms; by=x->x[2])[1:k])
	else
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		window = maximum([1, i-w]):minimum([length(X), i+w])
		norms[window] .= [(typemax(Int), Inf) for _ in window]
		return first.(sort(norms; by=x->x[2])[1:k])
	end
end

function nearest_neighbors_forecasting(x, R; k=1, w=5)
	return mean(R[nearest_neighbors(x, R; k=k, w=w) .+ 1])
end

u0 = zeros(2)
tr = Vector{Vector{Float64}}(undef, 10_000)
tr[1] = u0
for i in 2:10_000
	tr[i] = henon(tr[i-1])
end

h = first.(tr)

τ = 20 #plot(1:30, selfmutualinfo(h, 1:30))
d = 6  #plot(1:10, delay_afnn(h, τ, 1:10; w=5))

R = delay_embed(h, τ, d)

prediction = similar(R, 100)
for i in 1:100
	q = i == 1 ? R[end] : prediction[i-1]
	q̃ = nearest_neighbors_forecasting(q, R; k=3)
	prediction[i] = q̃
end
