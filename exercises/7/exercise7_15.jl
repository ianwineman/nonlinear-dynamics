using DifferentialEquations
using Plots
using Associations

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

plots = []

for k in [[0.075, 0.0], [0.05, 0.05], [0.0375, 0.0375]]
	println(k)
	prob = ODEProblem(coupledRössler!, u0, (0.0, 100.0), vcat(p0, k))
	solv = solve(prob)

	p1 = plot(
		solv, 
		idxs=(1),
		label="\$ x_1\$"
	)
	plot!(
		solv, 
		idxs=(4), 
		label="\$ y_1\$", 
		xlabel="\$ t\$",
		title="\$ (k_x, k_y) = ($(k[1]), $(k[2]))\$",
		titlefontsize=10,
		titleposition=:left
	)

	precise = true
	discretization = CodifyVariables(TransferOperator(RectangularBinning(5.0, precise)))
	est_disc_to = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
	te_xy = association(est_disc_to, getindex.(solv.u, 1), getindex.(solv.u, 4))
	te_yx = association(est_disc_to, getindex.(solv.u, 4), getindex.(solv.u, 1))

	te_xy = round(te_xy; digits=3)
	te_yx = round(te_yx; digits=3)

	p2 = plot(
		solv, 
		idxs=(1,4),
		label=false,
		xlabel="\$ x_1\$",
		ylabel="\$ y_1\$",
		legendfontsize=6,
		lw=0.5,
		aspect_ratio=:equal,
		title="\$ TE: x_1 → y_1 = $te_xy, y_1 → x_1 = $te_yx\$",
		titlefontsize=10,
		titleposition=:left
	)
	push!(plots, p1, p2)
end

plot(
	plots..., 
	layout=@layout[
		a{0.7w} b{0.3w}
		c{0.7w} d{0.3w}
		e{0.7w} f{0.3w}
	], 
	size=(900, 600), 
	plot_title="Coupled Rössler system", 
	plot_titlefontsize=12 
)
savefig("plots/exercise7_15.png")

