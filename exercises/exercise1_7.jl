using DifferentialEquations
using Plots

function single_pendulum!(du, u, p, t)
    θ, v = u
    du[1] = v
    du[2] = -p*sin(θ)
end

u0s = [[1.0,2.0], [1.0,3.0], [1.0,4.0]]
tspan = (0.0,10.0)
p = 1.0

sols = []
for u0=u0s
	prob = ODEProblem(single_pendulum!, u0, tspan, p)
	sol = solve(prob, dtmax=0.01)
	push!(sols,sol)
end

p1 = plot(
	[[sol.u[i][1] for i=1:length(sol.t)] for sol=sols],
	[[sol.u[i][2] for i=1:length(sol.t)] for sol=sols],
	xlabel="θ",
	ylabel="v",
	ylim=[0,5],
	title="Single pendulum",
	label=hcat(["(θ₀,v₀) = $(Tuple(sol.prob.u0))" for sol=sols]...),
	c=[:red :blue :orange],
	linewidth=2
)
p2 = plot(
	[sol.t for sol=sols],
	[[0.5*(sol.u[i][2])^2 - p*cos(sol.u[i][1]) for i=1:length(sol.t)] for sol=sols],
	xlabel="t",
	ylabel="H",
	ylim=[0,10],
	title="Hamiltonian of Single pendulum",
	label=hcat(["(θ₀,v₀) = $(Tuple(sol.prob.u0))" for sol=sols]...),
	c=[:red :blue :orange],
	linewidth=2
)

plot(p1,p2, layout=(2,1))
savefig("plots/exercise1_7.png")