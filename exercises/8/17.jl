using RecurrenceAnalysis
using Plots
using DifferentialEquations

function brusselator!(du, u, p, t)
   a, b = p
   x, y = u
   du[1] = 1 - (b+1)*x + a*(x^2)*y
   du[2] = b*x - a*(x^2)*y
end

colors = palette(:devon, 3) # purple, blue, white

prob = ODEProblem(brusselator!, [0.5909924611496271, 4.843385159662832], (0.0, 1_000.0), [0.3, 1.5])
solv = solve(prob, dtmax=0.1)

p1 = plot(
	solv, 
	idxs=(1, 2), 
	aspect_ratio=:equal,
	xlabel="\$ x\$",
	ylabel="\$ y\$",
	title="Brusselator",
	titlefontsize=10,
	label=false,
	lc=colors[2],
	ytick=2:6
)

len = 100
rt_averages = []
rt_entropies = []
dl_averages = []
dl_entropies = []
Ns, ϵs = range(100, 10_000; length=len), range(0.1, 0.05; length=len)
for (N, ϵ) in zip(Ns, ϵs)
	(N, ϵ) |> println
	X = solv.u[1:Int(N)]

	R = RecurrenceMatrix(StateSpaceSet(X), ϵ)
	push!(rt_averages, rt_average(R))
	push!(rt_entropies, rt_entropy(R))
	push!(dl_averages, dl_average(R))
	push!(dl_entropies, dl_entropy(R))
end

p2 = heatmap(
	RecurrenceMatrix(StateSpaceSet(solv.u[1:500]), 0.05),
	c=:devon,
	cbar=false,
	aspect_ratio=:equal,
	xlim=[1, 500],
	ylim=[1, 500],
	title="Recurrence matrix \$ (ϵ = 0.05)\$",
	titlefontsize=10,
)

p3 = plot(
	rt_averages ./ 10, # ./ 10 since dtmax=0.1,
	lc=colors[1],
	label="observed \$ ⟨r⟩\$",
	xlabel="\$ (N → ∞, ϵ → 0)\$",
	xlim=[1, 100]

)
hline!([12.558325799738558], lc=colors[2], label="limit cycle period") # 12.5... is limit cycle period
# rt_average -> period

p4 = plot(
	rt_entropies,
	xlabel="\$ (N → ∞, ϵ → 0)\$",
	ylabel="\$ H_r\$",
	label=false,
	lc=colors[2],
	xlim=[1, 100]
)
# rt_entropy -> 0 (all recurrences happen at period length) (converges slowly)

p5 = plot(
	dl_averages,
	xlabel="\$ (N → ∞, ϵ → 0)\$",
	label="observed \$ ⟨ℓ⟩\$",
	lc=colors[1],
	xlim=[1, 100]
)
plot!(1:len, log.(1:len) .+ 40, label="log", lc=colors[2])

p6 = plot(
	dl_entropies,
	xlabel="\$ (N → ∞, ϵ → 0)\$",
	ylabel="\$ H_ℓ\$",
	label=false,
	lc=colors[2],
	xlim=[1, 100]
)

plot(p1, p2, p3, p4, p5, p6, size=(600, 800), layout=(3, 2))
savefig("plots/8.17.png")

# dl_max -> ∞
# dt_average and dt_entropy should both grow without bound (roughly logarithmically)
# dt_average will grow since the number of arbitary length lines increases
# dt_entropy grows since the smallest line stays small and the distribution gets wider and wider
