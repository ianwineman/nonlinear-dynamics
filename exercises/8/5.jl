using Plots

lm(x) = 4 * x * (1 - x)
ρ(x) = 1 / (π * sqrt(x * (1 - x)))

X = [rand() for _ in 1:10_000]

for _ in 1:100
	global X = lm.(X)
end

colors = palette(:devon, 3) # purple, blue, white

histogram(X; 
	bins=25, 
	xlim=[0, 1], 
	ylim=[0, 3], 
	normalize=true, 
	label=false,
	xlabel="\$ x\$",
	title="Invariant density of Logistic map \$ (r=4)\$",
	c=colors[1]
)
plot!(
	range(0, 1; length=10_000), 
	ρ.(range(0, 1; length=10_000)), 
	label="\$ ρ(x)\$",
	lc=colors[2],
	lw=2
)
savefig("plots/8.5.png")
