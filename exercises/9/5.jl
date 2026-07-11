using DifferentialEquations
using Plots

function duffing!(du, u, p, t)
	a, d, ω = p

    du[1] = u[2]
    du[2] = -u[1] - (u[1])^3 - d*u[2] + a*sin(u[3])
    du[3] = ω
end

"return lorenz map of `x` for `x[m], x[m+k], x[m+2k], ...`"
function lorenzmap(x; m=1, k=1)
	[[x[m + i*k], x[m + (i+1)*k]] for i in 0:(((length(x) - m) ÷ k) - 1)]
end

T = 5.21
tmax = 100.0

prob = ODEProblem(duffing!, [1.0, 0.0, 0.0], (0.0, tmax), [19.75, 0.1, 2π/T])
solv = solve(prob)


p1 = plot(
	solv, idxs=(1, 2),
	label=false,
	xlabel="\$ x_1\$",
	ylabel="\$ x_2\$"
)

tmax = 10_000.0

prob = ODEProblem(duffing!, [1.0, 0.0, 0.0], (0.0, tmax), [19.75, 0.1, 2π/T])
solv = solve(prob)

pm = []
for n in 0:((tmax ÷ T) - 1 |> Int)
	t = T + n*T
	x1, x2, _ = solv(t)
	push!(pm, [x1, x2])
end

p2 = scatter(
	first.(pm), 
	last.(pm), 
	msw=0, 
	ms=1,
	label=false,
	xlabel="\$ x_1\$",
	ylabel="\$ x_2\$"
)

tmax = 1_000.0

prob = ODEProblem(duffing!, [1.0, 0.0, 0.0], (0.0, tmax), [19.75, 0.1, 2π/T])
solv = solve(prob)

pm = []
for n in 0:((tmax ÷ T) - 1 |> Int)
	t = T + n*T
	x1, x2, _ = solv(t)
	push!(pm, [x1, x2])
end

u1 = first.(pm)

p3 = plot(
	u1,
	label=false,
	xlabel="\$ t\$",
	ylabel="\$ x_1\$"
)

tmax = 10_000.0
prob = ODEProblem(duffing!, [1.0, 0.0, 0.0], (0.0, tmax), [19.75, 0.1, 2π/T])
solv = solve(prob, Vern9(), reltol=1e-9, abstol=1e-9, maxiters=10e7)

pm = []
for n in 0:((tmax ÷ T) - 1 |> Int)
	t = T + n*T
	x1, x2, _ = solv(t)
	push!(pm, [x1, x2])
end

u1 = first.(pm)

lm = lorenzmap(last.(u1); m=1501, k=4)

p4 = scatter(
	first.(lm), 
	last.(lm),
	label=false,
	msw=0,
	ms=1,
	xlabel="\$ x_1[n]\$",
	ylabel="\$ x_1[n+1]\$"
)

plot(p1, p2, p3, p4, plot_title="Duffing oscillator")
savefig("plots/9.5.png")
