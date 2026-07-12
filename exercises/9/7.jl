using DifferentialEquations
using Plots
using ProgressMeter

function toda!(du, u, p, t)
	a, d, ω = p

    du[1] = u[2]
    du[2] = a*sin(u[3]) - d*u[2] - exp(u[1]) + 1
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
        pm = poincarémap(2π/p[3], tmax, solv)
        if length(pm) > w
            period = unique(x -> round(x; digits=k), pm[(end-w):end] .|> first)
        else
            period = unique(x -> round(x; digits=k), pm .|> first)
        end
        for point in period
            push!(points, (p, point))
        end
    end
    return points
end

ω = 2.0
a = 0.88
T = 2π / ω
params = [a, 0.1, ω]
u0 = [1.0, 0.0, 0.0]
tmax = 100.0

prob = ODEProblem(toda!, u0, (0.0, tmax), params)
solv = solve(prob)

p1 = plot(
    solv, idxs=(1, 2),
    label="\$ ω=$ω\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$"
)

pm = poincarémap(T, tmax, solv)
p2 = scatter(
    first.(pm[(end-10):end]), 
    last.(pm[(end-10):end]),
    label="\$ ω=$ω\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$",
    ms=2,
    msw=0
)

p3 = plot(
    first.(pm[(end-10):end]),
    xlabel="\$ n\$",
    ylabel="\$ x_1[n]\$",
    label="\$ ω=$ω\$"
)

od = orbitdiagram(toda!, [1.0, 0.0, 0.0], [[0.88, 0.1, ω] for ω in 0.0:0.001:2.2]; tmax=tmax)
p4 = scatter(
    getindex.(first.(od), 3), 
    last.(od),
    msw=0,
    ms=0.5,
    label="\$ a=$a\$",
    xlabel="\$ ω\$",
    ylabel="\$ x_1\$"
)

lay = @layout [a b c; d]
plot(
    p1, p2, p3, p4,
    layout=lay, 
    plot_title="Toda oscillator", 
    plot_titlefontsize=10,
    size=(800, 400),
    margin=2.5Plots.mm
)
savefig("plots/9.7.png")
