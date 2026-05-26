using DynamicalBilliards
using Plots

# hyperparams:
widths    = 100   # length of width range
u0s_count = 1_000 # number of u0s per w
lst       = 100   # lyapunovspectrum t

chaotics = zeros(widths)
for (i, w) in enumerate(range(0.0, 2.0; length=widths))
	w == 0.0 && continue
	println(w)
	bm = billiard_mushroom(1.0, w, 1.0; door=false)
	u0s = [MushroomTools.randin_mushroom(1.0, w, 1.0) for _ in 1:u0s_count]

	chaotic = 0
	for u0 in u0s
		λ1 = lyapunovspectrum(u0, bm, lst)[1]
		if λ1 > 0.1
			chaotic += 1
		end
	end
	chaotics[i] = chaotic / u0s_count
end
plot(
	range(0.0, 2.0; length=widths), 
	chaotics,
	xlabel="\$ w\$",
	ylabel="Chaotic trajectory proportion",
	label=false,
	xlim=[0.0, 2.0],
	ylim=[0.0, 1.0],
	title="Mushroom billiard \$ (ℓ=1, w, r=1)\$"
)
savefig("plots/exercise8_9.png")
