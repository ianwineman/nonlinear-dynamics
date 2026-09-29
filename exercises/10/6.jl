using DifferentialEquations, Plots
using Plots.PlotMeasures
using Random

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

solvs = []
NRs = []
dumats = []
for (N, R) in zip(25:5:50, 0.35 .* collect(25:5:50) .|> ceil .|> Int)
	println("N=$N, R=$R")
	iseed = 5282 #7553,5 #3141,5 #1409,4 #9041,2
	rng = MersenneTwister(iseed)

	α = 1.46
	ω = 0.0
	params = [N, R, α, ω]
	u0 = 2 .* pi .* (rand(rng, N) .- 0.5)


	prob = ODEProblem(chimera!, u0, (0.0, 10000.0), params)
	solv = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=50.0)

	dus = map(zip(solv.u, solv.t)) do (u, t)
    du = similar(u)
	    chimera!(du, u, params, t)
	    du
	end

	dumat = reduce(hcat, dus)
	push!(solvs, solv)
	push!(NRs, (N, R))
	push!(dumats, dumat)
end

plots = []
for (i, s) in enumerate(solvs)
	N, R = NRs[i]
	push!(plots, heatmap(
		1:N, 
		s.t, 
		dumats[i]', 
		title="\$ N=$N, R=$R\$",
		ylabel="\$ t\$",
		xlabel="\$ i\$",
		formatter=:plain,
		xtick=[1, N]
	))
	ops = orderparameter(s)
	op_plot = plot(
		s.t, 
		ops, 
		title="\$ N=$N, R=$R\$", 
		label=false,
		xlabel="\$ t\$",
		ylabel="\$ Z(t)\$",
		ylim=[0, 1]
	)
	ff = findfirst(>(0.95), ops)
	if !isnothing(ff)
		scatter!(
			[s.t[ff]], 
			[ops[ff]], 
			label="\$ Z($(s.t[ff] |> Int)) > 0.95\$",
			msw=0,
			mc=:red
		)
	end
	push!(plots, op_plot)
end

plot(
	plots..., 
	size=(1200, 800),
	left_margin = 4mm,
    right_margin = 4mm,
    bottom_margin = 2mm,
    top_margin = 2mm
)
savefig("plots/10.6.png")
