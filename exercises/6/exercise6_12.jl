using Plots, Measures
using DynamicalSystems
using Statistics
using LinearAlgebra
using DifferentialEquations

function rössler!(du, u, p, t)
	a, b, c = p
	x, y, z = u

	du[1] = -y - z
	du[2] = x + a * y
	du[3] = b + z * (x - c)
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
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
	else
		norms = [(j, norm(x - y)) for (j, y) in enumerate(X)]
		window = maximum([1, i-w]):minimum([length(X), i+w])
		norms[window] .= [(typemax(Int), Inf) for _ in window]
		return [i for i in first.(sort(norms; by=x->x[2])[1:k]) if i < length(X)]
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

x_max_nrmse = []
for N in 200.0:100.0:5_000.0
	prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0, N), [0.2, 0.2, 5.7])
	solv = solve(prob; dtmax=0.1)

	h = first.(solv.u)

	τ = 12 # plot(1:30, selfmutualinfo(h, 1:30), xtick=1:30)
	d = 4  # plot(1:10, delay_afnn(h, τ, 1:10; w=5))

	R = delay_embed(h, τ, d)

	ground_truth = R[Int(N-99):end]
	prediction = similar(R, 100)
	for i in 1:100
		q = i == 1 ? R[Int(N - 100)] : prediction[i-1]
		q̃ = nearest_neighbors_forecasting(q, R; k=5) 
		push!(h, first(q̃))
		R = delay_embed(h, τ, d)
		prediction[i] = q̃
	end

	NRMSEs = zeros(100)
	for i in 1:100
		NRMSEs[i] = NRMSE(ground_truth[i], prediction[i], mean(ground_truth))
	end
	println("min/max NRMSE $((minimum(NRMSEs), maximum(NRMSEs))) [$N]")
	push!(x_max_nrmse, maximum(NRMSEs))
end

z_max_nrmse = []
for N in 200.0:100.0:5_000.0
	prob = ODEProblem(rössler!, [1, -2, 0.1], (0.0, N), [0.2, 0.2, 5.7])
	solv = solve(prob; dtmax=0.1)

	h = last.(solv.u)

	τ = 14 # plot(1:30, selfmutualinfo(h, 1:30), xtick=1:30)
	d = 5  # plot(1:10, delay_afnn(h, τ, 1:10; w=5))

	R = delay_embed(h, τ, d)

	ground_truth = R[Int(N-99):end]
	prediction = similar(R, 100)
	for i in 1:100
		q = i == 1 ? R[Int(N - 100)] : prediction[i-1]
		q̃ = nearest_neighbors_forecasting(q, R; k=5) 
		push!(h, first(q̃))
		R = delay_embed(h, τ, d)
		prediction[i] = q̃
	end

	NRMSEs = zeros(100)
	for i in 1:100
		NRMSEs[i] = NRMSE(ground_truth[i], prediction[i], mean(ground_truth))
	end
	println("min/max NRMSE $((minimum(NRMSEs), maximum(NRMSEs))) [N=$N]")
	push!(z_max_nrmse, maximum(NRMSEs))
end

plot(
	200.0:100.0:5_000.0,
    x_max_nrmse,
    xlabel="Time-series length \$N\$",
    ylabel="Maximum NRMSE (100 step prediction)",
    label="\$ h = x\$",
	title="Rössler system",
	lc=:blue,
	lw=1.5,
	ytick=vcat(1:2:13, 0.5)
)
plot!(
	200.0:100.0:5_000.0,
    z_max_nrmse,
    label="\$ h = z\$",
    lc=:orange,
    lw=1.5
)
hline!([0.5], ls=:dash, label="\$ n_c\$", lc=:black, lα=0.5)
savefig("plots/exercise6_12.png")
