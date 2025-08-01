using DifferentialEquations, Plots, FFTW, Statistics

function henonheiles!(du, u, p, t)
   x, y, vx, vy = u
   du[1] = vx
   du[2] = vy
   du[3] = -x - 2*x*y
   du[4] = -y - x^2 + y^2
end

u0s = (
	[0.0, -0.25, 0.42, 0.0],
	[0.0, 0.1, 0.5, 0.0],
	[0.0, 0.30266571044921875, 
	0.4205654433900762, 0.0],
)

plots = []
Δt = 0.05
for (i, u0) in enumerate(u0s)
	prob = ODEProblem(henonheiles!, u0, (0.0,1000.0))
	sol = solve(prob, dtmax=Δt)

	r = sol[1,:]
	P = abs2.(rfft(r .- mean(r)))
	ν = rfftfreq(length(r), 1/Δt)

	p = plot(ν, P ./ maximum(P))
	push!(plots, p)
end

plot(
	plots..., 
	layout=(3,1), 
	xlim=[0.0,0.4], 
	ylim=[10.0^(-5), 1.0], 
	yscale=:log10, 
	linewidth=3,
	linecolor=[:red :orange :blue],
	labels=["Chaotic" "Quasiperiodic" "Periodic"],
	xlabel=["" "" "Frequency, \$ν\$"],
	ylabel=["" "\$ P/max(P)\$" ""],
	plot_title="Hénon-Heiles normalized power spectra",
	plot_titlefontsize=14
)
savefig("plots/figure2_6.png")
