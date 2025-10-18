using Plots, DifferentialEquations

function brusselator!(du, u, p, t)
   a, b = p
   x, y = u
   du[1] = 1 - (b+1)*x + a*(x^2)*y
   du[2] = b*x - a*(x^2)*y
end

u0s = vcat(
   [[0.0,w] for w in 0.0:5.0],
   [[u,0.0] for u in 0.0:8.0],
   [[u,2.0*u + 5.0] for u in 0.0:0.5:1.0],
   [[u,-1.0*u + 8.0] for u in 1.0:8.0]
)
color_palette = palette(:winter, length(u0s))
plot(
   xlim=[-1,9], ylim=[-1,8], 
   title="Region \$ R\$ containing Brusselator limit cycle", 
   label=false,
   xlabel="\$ u\$",
   ylabel="\$ w\$"
)
for (i,u0) in enumerate(u0s)
   local prob = ODEProblem(brusselator!, u0, (0.0,40.0), [0.3, 1.5])
   local sol = solve(prob)
   plot!(sol, idxs=(1,2), label=false, linecolor=color_palette[i], linewidth=2)
end

plot!(0.0:0.1:8.0, [0.0 for _ in 0.0:0.1:8.0], linecolor=:grey50, linewidth=2, label="\$ δR\$")
plot!([0.0 for _ in 0.0:0.1:5.0], 0.0:0.1:5.0, linecolor=:grey50, linewidth=2, label=false)
plot!(0.0:0.1:1.0, 2 .* (0.0:0.1:1.0) .+ 5.0, linecolor=:grey50, linewidth=2, label=false)
plot!(1.0:0.1:8.0, -1 .* (1.0:0.1:8.0) .+ 8.0, linecolor=:grey50, linewidth=2, label=false)
savefig("plots/exercise2_2.png")
