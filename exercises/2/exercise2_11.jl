using Plots, Measures, DifferentialEquations; ENV["GKSwstype"] = "100"

function nose_hoover!(du, u, p, t)
   x, y, z = u
   du[1] = y
   du[2] = y*z - x
   du[3] = 1 - y^2
end

function trajectory(u0,tspan)
   prob = ODEProblem(nose_hoover!, u0, tspan)
   sol = solve(prob, dtmax=0.01)
   return sol
end

anim = @animate for t in 0.0:0.1:50.0
   sol1 = trajectory([0.0,1.5499334227822898,0.0],(0.0,t))
   sol2 = trajectory([0.0,1.0,0.0],(0.0,t))
   sol3 = trajectory([0.0,0.1,0.0],(0.0,t))

   p1 = plot(title="Nose-Hoover system", legendtitle="\$ t=$t\$")
   plot!(sol1, idxs=(1,2,3), linecolor=:red, label="Periodic", xlim=[-5,5], ylim=[-5,5], zlim=[-5,5])
   plot!(sol2, idxs=(1,2,3), linecolor=:orange, label="Quasiperiodic", xlim=[-5,5], ylim=[-5,5], zlim=[-5,5])
   plot!(sol3, idxs=(1,2,3), linecolor=:blue, label="Chaotic", xlim=[-5,5], ylim=[-5,5], zlim=[-5,5])
   plot(p1)


   println("t=$t")
end

gif(anim, "plots/exercise2_11.gif", fps=20);
