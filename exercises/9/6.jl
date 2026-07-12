using DifferentialEquations
using Plots
using ProgressMeter

function duffing!(du, u, p, t)
	a, d, ω = p

    du[1] = u[2]
    du[2] = -u[1] - (u[1])^3 - d*u[2] + a*sin(u[3])
    du[3] = ω
end

function poincarémap(T, tmax, solv)
    pm = []
    for n in 0:((tmax ÷ T) - 1 |> Int)
        t = T + n*T
        x1, x2, _ = solv(t)
        push!(pm, [x1, x2])
    end
    return pm
end

function orbitdiagram(system, u0, prange; i=1, tmax=100.0, w=10, k=1)
    points = []
    @showprogress for p in prange
        prob = ODEProblem(system, u0, (0.0, tmax), p)
        solv = solve(prob)
        pm = poincarémap(2π/p[3], tmax, solv)[(end-w):end]
        period = unique(x -> round(x; digits=k), pm .|> first)
        for point in period
            push!(points, (p, point))
        end
    end
    return points
end

#=
# Figure 9.2
od = orbitdiagram(duffing!, [0.1, 0.1, 0.0], [[7, 0.1, 2π/T] for T in 0.1:0.001:30]; tmax=1000.0, k=3)
scatter(
    2π ./ getindex.(first.(od), 3), 
    last.(od),
    msw=0,
    ms=1,
    mc=:blue,
    xtick=0:5:30
)
=#

T = 2π / 0.615
params = [12.0, 0.1, 0.615]
u0 = [1.0, 0.0, 0.0]
tmax = 1000.0

prob = ODEProblem(duffing!, u0, (0.0, tmax), params)
solv = solve(prob)

p1 = plot(
    solv, idxs=(1, 2),
    label="\$ ω=0.615\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$"
)

pm = poincarémap(T, tmax, solv)
p2 = scatter(
    first.(pm), 
    last.(pm),
    label="\$ ω=0.615\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$",
    ms=1,
    msw=0
)

od = orbitdiagram(duffing!, [1.0, 0.0, 0.0], [[12.0, 0.1, ω] for ω in 0.5:0.00001:0.7]; tmax=tmax, k=3)
p3 = vline([0.615], label=false)
scatter!(
    getindex.(first.(od), 3), 
    last.(od),
    msw=0,
    ms=0.5,
    label=false,
    xlabel="\$ ω\$",
    ylabel="\$ x_1\$",
    xtick=sort!(vcat(collect(0.5:0.05:0.7), 0.615))
)

lay = @layout [a b; c]
plot(
    p1, p2, p3, 
    layout=lay, 
    plot_title="Duffing oscillator [trajectory, Poincaré map; orbit diagram]", 
    plot_titlefontsize=10
)
savefig("plots/9.6.png")
