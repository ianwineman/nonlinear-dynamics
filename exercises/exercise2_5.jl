using Plots, Measures, DifferentialEquations

function lotka_volterra!(du, u, p, t)
   α, β, γ, δ = p
   x, y = u
   du[1] = α*x - β*x*y
   du[2] = γ*x*y - δ*y
end

α, β, γ, δ = 1.1, 0.4, 0.4, 0.1
sols = []

for y0 in [1.0, 2.0, 5.0, 7.0, 10.0, 12.0, 15.0]
   prob = ODEProblem(lotka_volterra!, [1.0,y0], (0.0,100.0), [α, β, γ, δ])
   sol = solve(prob, dtmax=0.01)
   push!(sols, sol)
end

p1 = plot(
   [[sol.u[i][1] for i=1:length(sol.t)] for sol in sols],
   [[sol.u[i][2] for i=1:length(sol.t)] for sol in sols],
   xlabel="x",
   ylabel="y",
   label=hcat(["$(sol.prob.u0[2])" for sol in sols]...),
   legendtitle="\$y_0\$",
)
vline!([δ/γ], linecolor=:black, linestyle=:dash, label="y-nullcline")
hline!([α/β], linecolor=:black, linestyle=:dashdot, label="x-nullcline")

p2 = plot(
   sols[3], 
   xlabel="time", 
   ylabel="population", 
   label=["x (prey)" "y (predator)"], 
   linecolor=[:blue :red],
   legendtitle="\$(x_0,y_0)=$(Tuple(sols[3].prob.u0))\$"
)
plot(p1, p2, size=(800,400), plot_title="Lotka-Volterra Equations\n\$ α, β, γ, δ = $α, $β, $γ, $δ\$", top_margin=10mm)
savefig("plots/exercise2_5.png")
