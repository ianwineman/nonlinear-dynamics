using DynamicalBilliards
using CairoMakie, Plots, LaTeXStrings

function random_particle()
	θ = rand() * 2π

	return Particle(SVector(5*cos(θ), 5*sin(θ)), SVector(-5*cos(θ), -5*sin(θ)))
end

function get_particle(θ)
	return Particle(SVector(5*cos(θ), 5*sin(θ)), SVector(-5*cos(θ), -5*sin(θ)))
end

function collision_count(x)
	length(unique(x[2:end]))
end

bd = Obstacle{Float64}[]

d1 = Disk(SVector(0.0, 1.0), 0.6)
d2 = Disk(SVector(cos(π/2 + (2π/3)), sin(π/2 + (2π/3))), 0.6)
d3 = Disk(SVector(cos(π/2 + 2*(2π/3)), sin(π/2 + 2*(2π/3))), 0.6)

push!(bd, d1, d2, d3)

billiard = Billiard(bd)

#fig, ax = bdplot(billiard)

#p1 = random_particle()
#x, y = DynamicalBilliards.timeseries!(p1, billiard, 20)
#lines!(ax, x, y; color=:red)

#fig
#save("plots/billiard_8.15.png", fig)


#=
particle_count = 1_000
collisions = zeros(Int, particle_count)
θs = []
for i in 1:particle_count
	x, y, vx, vy = DynamicalBilliards.timeseries!(random_particle(), billiard, 20)
	collisions[i] = collision_count(x)
	push!(θs, ((first(vx), first(vy)), (last(vx), last(vy))))
end

p1 = Plots.plot(
	1:particle_count, 
	collisions,
	label=false,
	xlabel="particle \$ i\$",
	ylabel="collision count",
	ylim=[1, maximum(collisions) + 1],
	ytick=1:(maximum(collisions) + 1)
)
=#

particle_count = 1_000
collisions = zeros(Int, particle_count)
Sϕo = zeros(particle_count)
for (i, ϕ) in enumerate(range(0, 2π/3; length=particle_count))
	x, y, vx, vy = DynamicalBilliards.timeseries!(get_particle(ϕ), billiard, 20)
	collisions[i] = collision_count(x)

	Sϕo[i] = acos(last(vx)) * sign(last(vy)) + π
end

p1 = Plots.plot(
	range(0, 2π/3; length=particle_count), 
	collisions,
	label=false,
	xlabel=L"\phi_i",
	ylabel="collision count",
	ylim=[1, maximum(collisions) + 1],
	ytick=1:(maximum(collisions) + 1)
)

p2 = Plots.scatter(
	range(0, 2π/3; length=particle_count), 
	Sϕo,
	label=false,
	msw=0,
	ms=2,
	xlim=[0, 2π/3],
	ylim=[0, 2π],
	xticks=([0, 2π/3], ["\$ 0\$", "\$ 2π/3\$"]),
	yticks=([0, 2π], ["\$ 0\$", "\$ 2π\$"]),
	xlabel=L"\phi_i",
	ylabel=L"\phi_o = \mathbf{S}(\phi_i)",
	leftmargin=5Plots.mm,
	rightmargin=5Plots.mm,
)

Plots.plot(p1, p2, plot_title="Chaotic scattering")
savefig("plots/8.15.png")
