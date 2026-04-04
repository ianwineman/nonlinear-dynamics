using DifferentialEquations
using Plots
using Statistics, LinearAlgebra
using DynamicalSystems

function delayembed(w::Vector{Float64}, τ::Int, d::Int)
	L = length(w) - (d-1)*τ
	masks = [[i + j*τ for j in 0:(d-1)] for i in 1:L]
	return [w[m] for m in masks]
end

function nearestneighbors(x, X; k=1, w=5)
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

function convergentcrossmapping(x, y, τ, d)
	mx = delayembed(x, τ, d)
	ψ = zeros(length(y))
	for (i, x) in enumerate(mx)
		nn = nearestneighbors(x, mx; k=5)
		u = [exp(-norm(x - mx[j]) / norm(x - mx[nn[end]])) for j in nn]
		w = u ./ sum(u)
		ψ[i] = sum(w .* y[nn])
	end
	return cor(y, ψ)
end

function coupledRössler!(du, u, p, t)
	a, b, c, wx, wy, kx, ky = p
	x1, x2, x3, y1, y2, y3  = u

    du[1] = -wx * x2 - x3
    du[2] = wx * x1 + a * x2 + kx * (y2 - x2)
    du[3] = b + x3 * (x1 - c)
    du[4] = -wy * y2 - y3
    du[5] = wy * y1 + a * y2 + ky * (x2 - y2)
    du[6] = b + y3 * (y1 - c)
end

u0 = [0.0, 0.2, 0.0, 0.11, 0.19, 0.1]
p0 = [0.2, 0.2, 5.7, 1.12, 1.08]

# x_1 for each (k_x, k_y)
xτ = [4, 4, 7] #plot(1:30, selfmutualinfo(getindex.(solv.u, 1), 1:30), xtick=1:30)
xd = [4, 6, 5] #plot(1:10, delay_afnn(getindex.(solv.u, 1), xτ[1], 1:10; w=5), xtick=1:10)

# y_1 for each (k_x, k_y)
yτ = [5, 5, 5] #plot(1:30, selfmutualinfo(getindex.(solv.u, 4), 1:30), xtick=1:30)
yd = [5, 4, 4] #plot(1:10, delay_afnn(getindex.(solv.u, 4), yτ[1], 1:10; w=5), xtick=1:10)

plots = []
for (i, k) in enumerate([[0.075, 0.0], [0.05, 0.05], [0.0375, 0.0375]])
	prob = ODEProblem(coupledRössler!, u0, (0.0, 100.0), vcat(p0, k))
	solv = solve(prob; dtmax=0.1)

	ccms = []
	for Ns in 50:50:length(solv.u)
		ccm = convergentcrossmapping(
			getindex.(solv.u, 1)[1:Ns], 
			getindex.(solv.u, 4)[1:Ns], 
			max(xτ[i], yτ[i]), 
			max(xd[i], yd[i])
		)
		println("$Ns $ccm")
		push!(ccms, (Ns, ccm))
	end
	println()
	p1 = plot(
		first.(ccms), 
		last.(ccms), 
		label="\$ x_1 → y_1\$",
		xlabel="\$ N_s\$",
		ylabel="\$ CCM\$",
		legendposition=:right,
		ylim=[0,1]
	)

	ccms = []
	for Ns in 50:50:length(solv.u)
		ccm = convergentcrossmapping(
			getindex.(solv.u, 4)[1:Ns], 
			getindex.(solv.u, 1)[1:Ns], 
			max(xτ[i], yτ[i]), 
			max(xd[i], yd[i])
		)
		println("$Ns $ccm")
		push!(ccms, (Ns, ccm))
	end
	println()
	plot!(
		first.(ccms), 
		last.(ccms), 
		label="\$ y_1 → x_1\$",
	)

	push!(plots, p1)
end

plot(
	plots..., 
	layout=(3, 1), 
	plot_title="Coupled Rössler system",
	plot_titlefontsize=12
)
savefig("plots/exercise7_16.png")
