using DifferentialEquations
using Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function plane(x, y; a=1, b=1, c=1, d=1)
	z = 30 - x
end


u0 = [2.0,1.0,1.0]
tspan = (0.0,500.0)
p = [10.0, 28.0, 8/3]
prob = ODEProblem(lorenz63!, u0, tspan, p)
sol = solve(prob)

p1 = plot(
	sol, idxs=(1,2,3),
	xlim=[-30,30], ylim=[-30,30], zlim=[0,60],
	xticks=[-30,0,30], yticks=[-30,0,30], zticks=[0,30,60],
	label="",
	linewidth=0.1,
	c=:red,
	xlabel="x", ylabel="y", zlabel="z",
	camera=[15,15]
)

x = -30:1:30
y = -30:1:30
X = [i for i in x, j in y]
Y = [j for i in x, j in y]
Z = 30 .- X

scatter!(
	X, Y, Z,
    label=:none,
	markerstrokewidth=0,
	markercolor=:grey,
	markeralpha=0.2,
	markersize=2
)

p2 = scatter(
	[u[1] for u=sol.u if isapprox(u[3],25.0,atol=0.1)],
	[u[2] for u=sol.u if isapprox(u[3],25.0,atol=0.1)],
	xlim=[-30,30], ylim=[-30,30],
	xticks=[-30,0,30], yticks=[-30,0,30],
	xlabel="x", ylabel="y",
	c=:red,
	markerstrokewidth=0,
	markersize=2,
	label="z=25.0",
	legendposition=:topleft
)

plot(
	p1, p2, layout=grid(1,2),
	plot_title="Lorenz-63 Poincaré Section"
)
savefig("plots/exercise1_9.png")