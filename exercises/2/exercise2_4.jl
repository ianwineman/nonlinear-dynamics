using DifferentialEquations
using Plots

function system!(du, u, p, t)
	x, y = u

	du[1] = y
	du[2] = -y * (x^2 + y^2 - 1) -x 
end

u0 = [0.001, 0.001]
tspan = (0.0,50.0)
p = []
prob = ODEProblem(system!, u0, tspan, p)
sol = solve(prob, dtmax=0.1)
plot(
	sol, idxs=(1,2), 
	aspect_ratio=:equal, 
	label="", xlabel="x", ylabel="y",
	color=:blue, lw=1.5,
	title="Limit Cycle in \$\\dot x = y, \\dot y = -y(x^2 + y^2 - 1) - x\$"
)

u0 = [-1.0, -1.0]
tspan = (0.0,50.0)
p = []
prob = ODEProblem(system!, u0, tspan, p)
sol = solve(prob, dtmax=0.1)
plot!(sol, idxs=(1,2), label="", color=:purple, lw=1.5)

u0 = [√2/2,√2/2]
tspan = (0.0,50.0)
p = []
prob = ODEProblem(system!, u0, tspan, p)
sol = solve(prob, dtmax=0.1)
plot!(sol, idxs=(1,2), lw=3, color=:black, label="")

savefig("plots/exercise2_4.png")