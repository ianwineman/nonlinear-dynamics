using Plots, Measures
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

function MSE(x, y)
	return sum((x .- y) .^ 2) / length(x)
end

function NRMSE(x, y, x̄)
	return sqrt(MSE(x, y) / MSE(x, x̄))
end

plots = []

ncs = []
for N in 10_000:5_000:100_000
	u0 = zeros(2)
	tr = Vector{Vector{Float64}}(undef, N + 100)
	tr[1] = u0
	for i in 2:(N + 100)
		tr[i] = henon(tr[i-1])
	end

	h = first.(tr)

	τ = 1 # discrete system
	d = 3 # plot(1:10, delay_afnn(h, τ, 1:10; w=5))

	R = delay_embed(h, τ, d)

	ground_truth = R[(N-99):end]
	prediction = similar(R, 100)
	for i in 1:100
		q = i == 1 ? R[N - 100] : prediction[i-1]
		q̃ = nearest_neighbors_forecasting(q, R; k=5) 
		push!(h, first(q̃))
		R = delay_embed(h, τ, d)
		prediction[i] = q̃
	end

	NRMSEs = zeros(100)
	for i in 1:100
		NRMSEs[i] = NRMSE(ground_truth[i], prediction[i], mean(ground_truth))
	end
	println("n_c = $(findfirst(>=(0.5), NRMSEs))")
	push!(ncs, findfirst(>=(0.5), NRMSEs))

	if N == 10_000
		p1 = plot(
			1:100, 
			NRMSEs, 
			ytick=0.0:0.5:5.0,
			lw=2,
			label=false,
			xlabel="Prediction length \$n\$",
			ylabel="\$ NRMSE\$",
			lc=:orange,
			xlim=[0,100],
			xtick=vcat(collect(0:25:100), [findfirst(>=(0.5), NRMSEs)]),
			left_margin=5mm,
			bottom_margin=5mm
		)
		hline!(
			[0.5], 
			label=false,
			ls=:dash,
			lc=:black,
			lα=0.25
		)
		scatter!(
			[findfirst(>=(0.5), NRMSEs)], 
			[0.5], 
			label="\$ n_c = $(findfirst(>=(0.5), NRMSEs))\$",
			mc=:red
		)
		push!(plots, p1)
	end
end
p2 = plot(
	10_000:5_000:100_000,
    ncs,
    xlabel="Time-series length \$N\$",
    ylabel="\$n_c\$",
    label=false,
    lw=2,
    lc=:green,
	left_margin=5mm,
	bottom_margin=5mm
)
scatter!([10_000], [ncs[1]], label="\$ n_c = $(ncs[1])\$", mc=:red)
push!(plots, p2)
plot(plots..., size=(800, 400), plot_title="Henon map")
savefig("plots/exercise6_11.png")
