using DynamicalBilliards
#using CairoMakie
#using Plots

function random_particle()
	ow = rand(1:4)

	if ow == 1     # top
		x = (rand() * 6) - 3
		y = 3.0

		θ = rand() * 2π
		u = cos(θ)
		v = -abs(sin(θ))
	elseif ow == 2 # right
		x = 3.0
		y = (rand() * 6) - 3

		θ = rand() * 2π
		u = -abs(cos(θ))
		v = sin(θ)
	elseif ow == 3 # bottom
		x = (rand() * 6) - 3
		y = -3.0

		θ = rand() * 2π
		u = cos(θ)
		v = abs(sin(θ))
	else           # left
		x = -3.0
		y = (rand() * 6) - 3

		θ = rand() * 2π
		u = abs(cos(θ))
		v = sin(θ)
	end

	u0 = SVector(x, y)
	ud = SVector(u, v)
	return Particle(u0, ud)
end

bd = Obstacle{Float64}[]

ow1 = FiniteWall(SVector(-3.0, 3.0), SVector(3.0, 3.0), SVector(0.0, -1.0))
ow2 = FiniteWall(SVector(-3.0, -3.0), SVector(3.0, -3.0), SVector(0.0, 1.0))
ow3 = FiniteWall(SVector(-3.0, -3.0), SVector(-3.0, 3.0), SVector(1.0, 0.0))
ow4 = FiniteWall(SVector(3.0, -3.0), SVector(3.0, 3.0), SVector(-1.0, 0.0))

iw1 = FiniteWall(SVector(-1.0, 1.0), SVector(1.0, 1.0), SVector(0.0, 1.0))
iw2 = FiniteWall(SVector(-1.0, -1.0), SVector(1.0, -1.0), SVector(0.0, -1.0))
iw3 = FiniteWall(SVector(-1.0, -1.0), SVector(-1.0, 1.0), SVector(-1.0, 0.0))
iw4 = FiniteWall(SVector(1.0, -1.0), SVector(1.0, 1.0), SVector(1.0, 0.0))

push!(bd, ow1, ow2, ow3, ow4, iw1, iw2, iw3, iw4)

billiard = Billiard(bd)

u0s = [random_particle() for _ in 1:1_000]
#Plots.scatter([first(p.pos) for p in u0s], [last(p.pos) for p in u0s])

#fig, ax = bdplot(billiard)

#x, y = DynamicalBilliards.timeseries!(u0s[1], billiard, 20)
#lines!(ax, x, y; color=:red)

#fig

nonchaotic = [DynamicalBilliards.lyapunovspectrum(u0, billiard, 100) == zeros(4) for u0 in u0s]
@assert nonchaotic == [1 for _ in 1:length(u0s)]
