using LinearAlgebra, DifferentialEquations, Plots

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

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

p = [10.0, 28.0, 8/3]
u0 = [20.0, 20.0, 20.0]

ϵcs = Vector{Tuple{Float64, Float64}}(undef, length(1.0:1.0:100.0))

Threads.@threads for (i,tend) in collect(enumerate(1.0:1.0:100.0))
	println(tend)
	prob = ODEProblem(lorenz63!, u0, (0.0,tend), p)
	solv = solve(prob)

	ϵs = [exp(j) for j in -20.0:0.1:5.0]
	Cs = [correlation_sum(solv.u, ϵ, 3) for ϵ in ϵs]
	ϵc = ϵs[findfirst(x -> !isinf(x), log.(Cs))]
	ϵcs[i] = (tend, ϵc)
end

plot(
	first.(ϵcs), 
	last.(ϵcs),
	xlim=[1,100],
	ylim=[0,maximum(last.(ϵcs))],
	lw=2,
	label=false,
	xlabel="\$t_{end}\$",
	ylabel="\$ϵ_c\$",
	title="Lorenz-63"
)
savefig("plots/exercise5_8.png")
