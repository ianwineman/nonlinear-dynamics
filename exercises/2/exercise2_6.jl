using Plots, Measures, DifferentialEquations

function lotka_volterra!(du, u, p, t)
   α, β, γ, δ = p
   x, y = u
   du[1] = α*x - β*x*y
   du[2] = γ*x*y - δ*y
end

function hamiltonian(x,y)
   return γ*x - δ*log(x) + β*y - α*log(y)
end

α, β, γ, δ = 1.1, 0.4, 0.1, 0.4

prob = ODEProblem(lotka_volterra!, [10.0,5.0], (0.0,50.0), [α, β, γ, δ])
sol = solve(prob, dtmax=0.01)

p1 = plot(
   sol,
   idxs=(1,2),
   xlabel="x",
   ylabel="y",
   label=false,
   legendtitle="\$(x_0,y_0)=$(Tuple(sol.prob.u0))\$",
   linecolor=:orange,
   linewidth=2,
   left_margin=5mm,
   right_margin=5mm
)
hline!([α/β], linecolor=:black, linestyle=:dashdot, label="x-nullcline")
vline!([δ/γ], linecolor=:black, linestyle=:dash, label="y-nullcline")

p2 = plot(
   sol, 
   xlabel="time", 
   ylabel="population", 
   label=["x (prey)" "y (predator)"], 
   linecolor=[:blue :red],
   linewidth=2,
   legendtitle="\$(x_0,y_0)=$(Tuple(sol.prob.u0))\$",
   legendposition=:topright,
   right_margin=5mm
)

p3 = plot(
   sol.t,
   [hamiltonian(u[1],u[2]) for u in sol.u],
   label="\$ H(x,y,t)\$",
   ylim=[0,1],
   linecolor=:black,
   linewidth=2,
   xlabel="time",
   ylabel="\$ H\$"
)
plot(p1, p2, p3, layout=(1,3), size=(1200,400), plot_title="Lotka-Volterra Equations\n\$ α, β, γ, δ = $α, $β, $γ, $δ\$", top_margin=10mm, bottom_margin=10mm)
savefig("plots/exercise2_6.png")
