using DifferentialEquations
using Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

u0 = [2.0,1.0,1.0]
tspan = (0.0,500.0)
p = [10.0, 28.0, 8/3]
prob = ODEProblem(lorenz63!, u0, tspan, p)
sol = solve(prob)

anim = @animate for z=0:0.1:50
	p1 = plot(
		sol, idxs=(1,2,3),
		xlim=[-30,30], ylim=[-30,30], zlim=[0,60],
		label="",
		linewidth=0.1,
		c=:orange,
		xlabel="x", ylabel="y", zlabel="z"
	)

	x = -30:0.1:30
	y = -30:0.1:30
	surface!(
		x, y, (x, y) -> z,
		alpha=0.3,
		c=:grey,
		colorbar=false
	)

	p2 = scatter(
		[u[1] for u=sol.u if isapprox(u[3],z,atol=0.1)],
		[u[2] for u=sol.u if isapprox(u[3],z,atol=0.1)],
		xlim=[-30,30], ylim=[-30,30],
		xlabel="x", ylabel="y",
		c=:orange,
		markerstrokewidth=0,
		markersize=2,
		label="z=$z",
		legendposition=:topleft
	)

	plot(
		p1, p2,
		plot_title="Lorenz-63"
	)

    println("z=$z")
end

gif(anim, "plots/exercise1_8.gif", fps=10)