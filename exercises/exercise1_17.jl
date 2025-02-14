using Plots
using DifferentialEquations

function SIR!(du, u, p, t)
	β, γ = p(t)
	S, I, R = u

	du[1] = -β*((S*I)/(S+I+R))
	du[2] = β*((S*I)/(S+I+R)) - γ*I
	du[3] = γ*I
end

u0_1 = [999.0, 1.0, 0.0]
u0_2 = [900.0, 100.0, 0.0]
tspan = (0.0,50.0)
p = [0.5, 0.05]

prob_1 = ODEProblem(SIR!, u0_1, tspan, t -> (0.5, 0.05))
sol_1 = solve(prob_1, dtmax=0.01)

prob_2 = ODEProblem(SIR!, u0_2, tspan, t -> (0.5, 0.05))
sol_2 = solve(prob_2, dtmax=0.01)

prob_3 = ODEProblem(SIR!, u0_1, tspan, t -> (10<t<30) ? (0.1, 0.05) : (0.5, 0.05))
sol_3 = solve(prob_3, dtmax=0.01)

prob_4 = ODEProblem(SIR!, u0_2, tspan, t -> (10<t<30) ? (0.1, 0.05) : (0.5, 0.05))
sol_4 = solve(prob_4, dtmax=0.01)

p1 = plot(
	sol_1.t, 
	[[u[1] for u=sol_1.u], 
     [u[2] for u=sol_1.u], 
     [u[3] for u=sol_1.u],
     [u[1]+u[2]+u[3] for u=sol_1.u]],
    label=["Susceptible" "Infected" "Recovered" "Total"],
    legendposition=:topright,
    xlabel="Time, t", ylabel="Count",
    xlim=[0,50], ylim=[0,1020],
    title="\$(S_0,I_0,R_0)=($(Int(sol_1.prob.u0[1])),$(Int(sol_1.prob.u0[2])),$(Int(sol_1.prob.u0[3])))\$, \$β=$(sol_1.prob.p(0)[1]), γ=$(sol_1.prob.p(0)[2])\$",
    titlefontsize=10,
    linewidth=[1.5 1.5 1.5 1.5], linecolor=[:gold :red2 :green3 :grey50]
)

p2 = plot(
	sol_2.t, 
	[[u[1] for u=sol_2.u], 
     [u[2] for u=sol_2.u], 
     [u[3] for u=sol_2.u],
     [u[1]+u[2]+u[3] for u=sol_2.u]],
    label=["Susceptible" "Infected" "Recovered" "Total"],
    legendposition=:topright,
    xlabel="Time, t", ylabel="Count",
    xlim=[0,50], ylim=[0,1020],
    title="\$(S_0,I_0,R_0)=($(Int(sol_2.prob.u0[1])),$(Int(sol_2.prob.u0[2])),$(Int(sol_2.prob.u0[3])))\$, \$β=$(sol_2.prob.p(0)[1]), γ=$(sol_2.prob.p(0)[2])\$",
    titlefontsize=10,
    linewidth=[1.5 1.5 1.5 1.5], linecolor=[:gold :red2 :green3 :grey50],
)

p3 = plot(
	sol_1.t, 
	[[u[1] for u=sol_3.u], 
     [u[2] for u=sol_3.u], 
     [u[3] for u=sol_3.u],
     [u[1]+u[2]+u[3] for u=sol_3.u]],
    label=["Susceptible" "Infected" "Recovered" "Total"],
    legendposition=:topright,
    xlabel="Time, t", ylabel="Count",
    xlim=[0,50], ylim=[0,1020],
    title="Quarantine for \$t∈[10,30], β=0.1\$",
    titlefontsize=10,
    linewidth=[1.5 1.5 1.5 1.5], linecolor=[:gold :red2 :green3 :grey50]
)

p4 = plot(
	sol_2.t, 
	[[u[1] for u=sol_4.u], 
     [u[2] for u=sol_4.u], 
     [u[3] for u=sol_4.u],
     [u[1]+u[2]+u[3] for u=sol_4.u]],
    label=["Susceptible" "Infected" "Recovered" "Total"],
    legendposition=:topright,
    xlabel="Time, t", ylabel="Count",
    xlim=[0,50], ylim=[0,1020],
    title="Quarantine for \$t∈[10,30], β=0.1\$",
    titlefontsize=10,
    linewidth=[1.5 1.5 1.5 1.5], linecolor=[:gold :red2 :green3 :grey50],
)

plot(
	p1, p2, p3, p4, 
	size=(800,800), 
	layout=(2,2),
	plot_title="SIR Model",
	plot_titlefontsize=12
)

savefig("plots/exercise1_17.png")