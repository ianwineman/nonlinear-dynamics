using DifferentialEquations, Plots, PlotlyJS
plotlyjs()

function gissinger!(du, u, p, t)
    μ, ν, γ = p
    x, y, z = u

    du[1] = μ*x - y*z
    du[2] = -ν*y + x*z
    du[3] = γ - z + x*y
end

sols = []
for μ=0.109:0.01:0.129
    u0 = [3.0,3.0,3.0]
    tspan = (0.0,250.0)
    p = [μ,0.1,0.9]
    prob = ODEProblem(gissinger!, u0, tspan, p)
    sol = solve(prob)
    push!(sols,sol)
end

p = Plots.plot(
    xlim=(-10,10), ylim=(-10,10), zlim=(-10,10),
    xguide="X", yguide="Y", zguide="Z",
    aspect_ratio=:equal 
)

for sol=sols
    plot!(
        [u[1] for u in sol.u], 
        [u[2] for u in sol.u], 
        [u[3] for u in sol.u], 
        label="μ=$(sol.prob.p[1])"
    )
end

display(p)