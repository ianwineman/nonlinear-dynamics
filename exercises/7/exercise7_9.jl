using Statistics, LinearAlgebra
using DifferentialEquations
using DynamicalSystems
using Plots

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

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end

prob = ODEProblem(lorenz63!, [20.0, 20.0, 20.0], (0.0,200.0), [10.0, 40.0, 8/3])
solv = solve(prob)

τ = 8 #plot(1:30, selfmutualinfo(first.(solv.u), 1:30), xtick=1:30)
d = 6 #plot(1:10, delay_afnn(first.(solv.u), τ, 1:10; w=5))

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 1), getindex.(solv.u[1:Ns], 2), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot(
	first.(ccms), 
	last.(ccms), 
	label="\$ x → y\$",
	xlabel="\$ N_s\$",
	ylabel="CCM measure",
	title="Lorenz 63",
	legendposition=:right
)

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 2), getindex.(solv.u[1:Ns], 1), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(first.(ccms), last.(ccms), label="\$ y → x\$")

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 1), getindex.(solv.u[1:Ns], 3), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(first.(ccms), last.(ccms), label="\$ x → z\$")

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 3), getindex.(solv.u[1:Ns], 1), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(first.(ccms), last.(ccms), label="\$ z → x\$")

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 2), getindex.(solv.u[1:Ns], 3), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(first.(ccms), last.(ccms), label="\$ y → z\$")

ccms = []
for Ns in 30:100:3000
	ccm = convergentcrossmapping(getindex.(solv.u[1:Ns], 3), getindex.(solv.u[1:Ns], 2), τ, d)
	println("$Ns $ccm")
	push!(ccms, (Ns, ccm))
end
plot!(first.(ccms), last.(ccms), label="\$ z → y\$")
savefig("plots/exercise7_9.png")
