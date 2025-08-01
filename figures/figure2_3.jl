using Plots, DifferentialEquations

function brusselator!(du, u, p, t)
   a, b = p
   x, y = u
   du[1] = 1 - (b+1)*x + a*(x^2)*y
   du[2] = b*x - a*(x^2)*y
end

prob = ODEProblem(brusselator!, [0.5909924611496271,4.843385159662832], (0.0,20.0), [0.3, 1.5])
sol = solve(prob)

limit_cycle_period = 12.558325799738558
lc_us = [
	solve(
		ODEProblem(brusselator!, [0.5909924611496271,4.843385159662832], (0.0,t), [0.3, 1.5])
	)[end][1] 
	for t=0.0:(limit_cycle_period/10):limit_cycle_period
]
lc_ws = [
	solve(
		ODEProblem(brusselator!, [0.5909924611496271,4.843385159662832], (0.0,t), [0.3, 1.5])
	)[end][2] 
	for t=0.0:(limit_cycle_period/10):limit_cycle_period
]

plot(
	sol, idxs=(1,2),
	xlabel="\$u\$", ylabel="\$w\$", label="",
	xlim=[0.5,2.2], ylim=[2.5,6.0],
	xtick=[1,2], ytick=[3,5],
	lw=2, lc=:slategray1,
	title="Brusselator, \$a=0.3,b=1.5\$",
)
scatter!(
	[1.0], [5.0],
	label="Unstable fixed point: \$(1,b/a)\$",
	mc=:white
)
scatter!(
	lc_us, lc_ws,
	label="",
	marker_z=[(2π*t)/limit_cycle_period for t=0.0:(limit_cycle_period/10):limit_cycle_period],
	mc=:reds,
	cbar_title="\$ϕ\$"
)

savefig("plots/figure2_3.png")