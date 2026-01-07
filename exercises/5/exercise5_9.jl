using DifferentialEquations, LinearAlgebra
using Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function neg_to_pos(vec)
	indices = []
	neg = findall(x->x<0, vec)
	pos = findall(x->x>0, vec)

	for i in neg
		if (i + 1) in pos
			push!(indices, (i, i+1))
		end
	end
	return indices
end

function psos(system, u0, p, tspan; plane=1, transient_ratio=0.75)
	prob = ODEProblem(system, u0, tspan, p)
	solu = solve(prob, dtmax=0.01)

	cross_plane = neg_to_pos(getindex.(solu.u, plane))
	ts = [solu.t[pair[1]] + (solu.t[pair[2]]-solu.t[pair[1]])/2 for pair in cross_plane]
	return [solu(t) for t in filter(x->x>tspan[2]*transient_ratio, ts)]
end

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

p  = [10.0, 28.0,  8/3]
u0 = [20.0, 20.0, 20.0]
tspan = (0.0, 200.0)

ps = psos(lorenz63!, u0, p, tspan; transient_ratio=0.25)

scatter(getindex.(ps, 1), getindex.(ps, 2), getindex.(ps, 3))

ϵs = [exp(i) for i in -20.0:0.1:5.0]
Cs = [correlation_sum(ps, ϵ, 3) for ϵ in ϵs]
n = findfirst(x->!isinf(x) && !isnan(x), log.(Cs))
plot(
	log.(ϵs)[n:end], 
	log.(Cs)[n:end],
	xlabel="\$\\log{(ϵ)}\$",
	xtick=-4:1:5,
	ylabel="\$\\log{(C)}\$",
	title="Lorenz-63 Poincaré Surface",
	label=false,
	lw=2,
	lc=:blue
)

xstart, xend = -1.0, 2.0
dim = (log.(Cs)[findfirst(x->x==xend, -20.0:0.1:5.0)] - log.(Cs)[findfirst(x->x==xstart, -20.0:0.1:5.0)])/(xend - xstart)
plot!(
	range(xstart,xend; length=100),
	range(xstart,xend; length=100) .* dim .+ log.(Cs)[findfirst(x->x==0.0, -20.0:0.1:5.0)] .- 0.3,
	label="\$Δ^{(C)} ≈ $(round(dim; digits=2))\$",
	lc=:lightcoral,
	lw=2
)
savefig("plots/exercise5_9.png")
