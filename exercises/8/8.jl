using DynamicalBilliards
using CairoMakie, Plots

bd = Obstacle{Float64}[]
sc1 = Semicircle(SVector(0.0, 0.5), 0.5, SVector(1.0, 0.0))
sc2 = Semicircle(SVector(1.0, 0.5), 0.5, SVector(-1.0, 0.0))

w1 = FiniteWall(SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0))

δ = 0.025 # stem width
w2 = FiniteWall(SVector(0.0, 1.0), SVector(0.5 - δ, 1.0), SVector(0.0, -1.0))
w3 = FiniteWall(SVector(1.0, 1.0), SVector(0.5 + δ, 1.0), SVector(0.0, -1.0))

w4 = FiniteWall(SVector(0.5 - δ, 1.0), SVector(0.5 - δ, 1.2), SVector(1.0, 0.0))
w5 = FiniteWall(SVector(0.5 + δ, 1.0), SVector(0.5 + δ, 1.2), SVector(-1.0, 0.0))

w6 = FiniteWall(SVector(0.5 - δ, 1.2), SVector(-0.2, 1.2), SVector(0.0, 1.0))
w7 = FiniteWall(SVector(0.5 + δ, 1.2), SVector(1.2, 1.2), SVector(0.0, 1.0))

w8 = FiniteWall(SVector(-0.2, 1.2), SVector(-0.2, 2.7), SVector(1.0, 0.0))
w9 = FiniteWall(SVector(1.2, 1.2), SVector(1.2, 2.7), SVector(-1.0, 0.0))

w10 = FiniteWall(SVector(-0.2, 2.7), SVector(1.2, 2.7), SVector(0.0, -1.0))

push!(bd, sc1, sc2, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10)

billiard = Billiard(bd)

u0s = [Particle(SVector(0.5, 0.5), SVector(cos(rand() * 2π), sin(rand() * 2π))) for _ in 1:10_000]

#fig, ax = bdplot(billiard)

#x, y = DynamicalBilliards.timeseries!(u0s[1], billiard, 200)
#lines!(ax, x, y; color=:blue)

#save("plots/billiard_8.8.png", fig)

sim_length = 1_000
above_one_count = zeros(sim_length + 1)
for u0 in u0s
	x, y = DynamicalBilliards.timeseries!(u0s[1], billiard, sim_length)

	for i in 1:(sim_length + 1)
		if y[i] > 1.0
			above_one_count[i] += 1
		end
	end
end

Plots.plot(
	1:(sim_length + 1), 
	above_one_count ./ 10_000,
	xlabel="collisions",
	ylabel="proportion of balls above \$ y=1\$",
	xlim=[0, sim_length],
	leftmargin=5Plots.mm,
	rightmargin=5Plots.mm,
	label=false
)
savefig("plots/8.8.png")
