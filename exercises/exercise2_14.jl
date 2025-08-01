using DifferentialEquations, Plots, FFTW, Statistics

function henonheiles!(du, u, p, t)
   x, y, vx, vy = u
   du[1] = vx
   du[2] = vy
   du[3] = -x - 2*x*y
   du[4] = -y - x^2 + y^2
end

u0s = (
	[0.0, -0.31, 0.354, 0.059],
	[0.0, 0.1, 0.5, 0.0],
	[0.0, -0.091, 0.46, -0.173],
)

plots = []
Δt = 0.05
for (i, u0) in enumerate(u0s)
	prob = ODEProblem(henonheiles!, u0, (0.0,1000.0))
	sol = solve(prob, dtmax=Δt)

	r = sol[1,:]
	P = abs2.(rfft(r .- mean(r)))
	ν = rfftfreq(length(r), 1/Δt)

	p = plot(ν, P ./ maximum(P), label=false, legendposition=:topright)
	push!(plots, p)
end

vline!(
	plots[1],
	[0.035],
	label="\$ v_1 ≈ 0.035\$"
)
vline!(
	plots[1],
	[0.123],
	label="\$ v_2 ≈ 0.123\$"
)
vline!(
	plots[1],
	[0.158],
	label="\$ v_1 + v_2\$"
)
vline!(
	plots[1],
	[0.245],
	label="\$ 2v_2\$"
)
vline!(
	plots[1],
	[0.315],
	label="\$ 2v_1 + 2v_2\$"
)
vline!(
	plots[2],
	[0.047],
	label="\$ v_1 ≈ 0.047\$"
)
vline!(
	plots[2],
	[0.114],
	label="\$ v_2 ≈ 0.114\$"
)
vline!(
	plots[2],
	[0.161],
	label="\$ v_1 + v_2\$"
)
vline!(
	plots[2],
	[0.229],
	label="\$ 2v_2\$"
)
vline!(
	plots[2],
	[0.322],
	label="\$ 2v_1 + 2v_2\$"
)
vline!(
	plots[3],
	[0.034],
	label="\$ v_1 ≈ 0.034\$"
)
vline!(
	plots[3],
	[0.137],
	label="\$ v_2 ≈ 0.137\$"
)
vline!(
	plots[3],
	[0.171],
	label="\$ v_1 + v_2\$"
)
vline!(
	plots[3],
	[0.274],
	label="\$ 2v_2\$"
)
vline!(
	plots[3],
	[0.308],
	label="\$ v_1 + 2v_2\$"
)

plot(
	plots..., 
	layout=(3,1), 
	xlim=[0.0,0.45], 
	xticks=collect(0.0:0.05:0.45),
	minorgrid=true,
	ylim=[10.0^(-5), 1.0], 
	yscale=:log10, 
	legendtitle=["\$ W≈3.514\$" "\$ W≈2.426\$" "\$ W≈4.029\$"],
	linestyle=[:solid :dot :dot :dot :dot :dot :solid :dot :dot :dot :dot :dot :solid :dot :dot :dot :dot :dot],
	linecolor=[:grey :red :blue :green3 :darkorange :purple3 :grey :red :blue :green3 :darkorange :purple3 :grey :red :blue :green3 :darkorange :purple3],
	xlabel=["" "" "Frequency, \$ν\$"],
	ylabel=["" "\$ P/max(P)\$" ""],
	plot_title="Hénon-Heiles normalized power spectra",
	plot_titlefontsize=14,
	size=(900,600)
)
savefig("plots/exercise2_14.png")
