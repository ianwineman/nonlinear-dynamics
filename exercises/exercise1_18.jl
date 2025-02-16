using DifferentialEquations
using Plots

function system!(du, u, p, t)
    v, θ = u
    du[1] = -sin(θ)
    du[2] = (-cos(θ)+v^2)/v
end

u0s = [[1.0,2.0], [1.0,3.0], [1.0,4.0]]
tspan = (0.0,25.0)

sols = []
for u0=u0s
	local prob = ODEProblem(system!, u0, tspan, p)
	local sol = solve(prob, dtmax=0.01)
	push!(sols,sol)
end

p1 = plot(
	[[sol.u[i][1] for i=1:length(sol.t)] for sol=sols],
	[[sol.u[i][2] for i=1:length(sol.t)] for sol=sols],
	xlabel="v",
	ylabel="θ",
	xlim=[0,3],
	ylim=[0,40],
	title="System: \$\\dot v = \\sin θ, \\dot θ=(v^2 - \\cos θ)/v\$",
	label=hcat(["(v₀,θ₀) = $(Tuple(sol.prob.u0))" for sol=sols]...),
	c=[:red :blue :orange],
	linewidth=2,
	legend=:bottomright
)
p2 = plot(
	[sol.t for sol=sols],
	[[sol.u[i][1] * cos(sol.u[i][2]) - (1/3)*sol.u[i][1]^3 for i=1:length(sol.t)] for sol=sols],
	xlabel="t",
	ylabel="H",
	ylim=[-2,0],
	title="Hamiltonian: \$H = v\\cos θ-\\frac{1}{3}v^3\$",
	label=hcat(["(v₀,θ₀) = $(Tuple(sol.prob.u0))" for sol=sols]...),
	c=[:red :blue :orange],
	linewidth=2,
	legend=:bottomright
)

plot(p1,p2, layout=grid(2,1; heights=(0.7,0.3)), size=(400,600))

savefig("plots/exercise1_18.png")