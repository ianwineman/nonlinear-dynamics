using DifferentialEquations
using Plots
using ProgressMeter

function vanderpol!(du, u, p, t)
	a, d, ω = p

    du[1] = u[2]
    du[2] = a*sin(u[3]) - d*((u[1])^2 - 1)*u[2] - u[1]
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
        solv = solve(prob, Vern9(), reltol=1e-9, abstol=1e-9, maxiters=10e7)
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

ω = 2.45
T = 2π / ω
params = [5.0, 5.0, ω]
u0 = [1.0, 0.0, 0.0]
tmax = 100.0

prob = ODEProblem(vanderpol!, u0, (0.0, tmax), params)
solv = solve(prob, Vern9(), reltol=1e-9, abstol=1e-9, maxiters=10e7)

p1 = plot(
    solv, idxs=(1, 2),
    label="\$ ω=$ω\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$"
)

pm = poincarémap(T, tmax, solv)
p2 = scatter(
    first.(pm[10:end]), 
    last.(pm[10:end]),
    label="\$ ω=$ω\$",
    xlabel="\$ x_1\$",
    ylabel="\$ x_2\$",
    ms=2,
    msw=0
)

p3 = plot(
    first.(pm[10:end]),
    xlabel="\$ n\$",
    ylabel="\$ x_1[n]\$",
    label="\$ ω=$ω\$"
)

small_range = 2.45:0.000001:2.475
od = orbitdiagram(vanderpol!, [1.0, 0.0, 0.0], [[5.0, 5.0, ω] for ω in 2.0:0.0001:5.0]; tmax=tmax)

od_wide = [o for o in od if o[1][3] ∉ small_range]
od_small_range = [o for o in od if o[1][3] ∈ small_range]

p4 = scatter(
    getindex.(first.(od_wide), 3), 
    last.(od_wide),
    msw=0,
    ms=0.5,
    label=false,
    xlabel="\$ ω\$",
    ylabel="\$ x_1\$"
)
scatter!(
    getindex.(first.(od_small_range), 3), 
    last.(od_small_range),
    msw=0,
    ms=0.5,
    mc=:orange,
    label=false
)

od = orbitdiagram(vanderpol!, [1.0, 0.0, 0.0], [[5.0, 5.0, ω] for ω in small_range]; tmax=tmax)
p5 = scatter(
    getindex.(first.(od), 3), 
    last.(od),
    msw=0,
    ms=0.5,
    label=false,
    xlabel="\$ ω\$",
    ylabel="\$ x_1\$",
    mc=:orange
)

lay = @layout [a b c; d; e]
plot(
    p1, p2, p3, p4, p5,
    layout=lay, 
    plot_title="van der Pol oscillator", 
    plot_titlefontsize=10,
    size=(800, 600),
    margin=2.5Plots.mm
)
savefig("plots/9.8.png")
