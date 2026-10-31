using DifferentialEquations, Plots
using Random
using LaTeXStrings
using Statistics

function chimera!(du, u, p, t)
	N, R, α, ω = p

	N, R = Int.([N, R])
    for i in 1:N
    	du[i] = ω - (1/(2*R)) * sum([sin(u[i] - u[mod(j, 1:N)] + α) for j in i-R:i+R])
    end
end

function orderparameter(solv, t)
	u = solv(t)
	abs(sum(exp.(u .* im)) / length(u))
end

function orderparameter(solv, ts::StepRange)
	[orderparameter(solv, t) for t in ts]
end

function orderparameter(solv)
	[abs(sum(exp.(u .* im)) / length(u)) for u in solv.u]
end


iseed = 5282 #7553,5 #3141,5 #1409,4 #9041,2
rng = MersenneTwister(iseed)

N = 40
R = 14
#α = 1.46
ω = 0.0
#params = [N, R, α, ω]
u0 = 2 .* pi .* (rand(rng, N) .- 0.5)

αs = sort(vcat(range(1.4, 1.52; length=50), 1.46))

maxZs = zeros(length(αs))

Threads.@threads for (i, α) in collect(enumerate(αs))
	@show i
	params = [N, R, α, ω]
	prob = ODEProblem(chimera!, u0, (0.0, 10000.0), params)
	solv = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=50.0)
	maxZs[i] = maximum(orderparameter(solv))
end

hspan([0.95, 1.0], color=:orange, alpha=0.1, label="sync"); vline!([1.46], lc=:orange, label=false); plot!(
	αs, 
	maxZs, 
	marker=:circ,
	xlabel="\$ α\$",
	ylabel=L"\max_{t} \;\; Z(t)",
	label=false,
	lc=:blue,
	mc=:blue,
	ms=2,
	xtick=1.4:0.01:1.52
)
savefig("plots/10.7.png")
