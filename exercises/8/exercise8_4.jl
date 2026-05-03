using LinearAlgebra
using Plots
using Statistics

struct Billiard
	xrange::UnitRange{Int64}
	yrange::UnitRange{Int64}
	circr::Float64
end

function nextcollision(
	B::Billiard, 
	δBxy::Tuple{Float64, Float64}, 
	unitdir::Tuple{Float64, Float64}
)
	x, y = δBxy
	u, v = unitdir
	top = last(B.yrange)
	bot = first(B.yrange)
	left = first(B.xrange)
	right = last(B.xrange)
	r = B.circr

	ttop = (top - y) / v   # collision time to y = top
	tbot = (bot - y) / v   # collision time to y = bot
	tlef = (left - x) / u  # collision time to x = left
	trig = (right - x) / u # collision time to x = right

	coltimes = [ttop, tbot, tlef, trig]
	t = minimum(filter(>(0), coltimes))
	δi = findfirst(==(t), coltimes)

	xx = x + (t * u)
	yy = y + (t * v)
	if δi == 1
		m = (top - y) / (xx - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with top
			return t, (xx, top), (u, -v), false
		else     # collision with circle
			x1  = (-m*b + √Δ) / (m^2 + 1)
			x2  = (-m*b - √Δ) / (m^2 + 1)
			y11 =  sqrt(r^2 - x1^2)
			y12 = -sqrt(r^2 - x1^2)
			y21 =  sqrt(r^2 - x2^2)
			y22 = -sqrt(r^2 - x2^2)
			points = [
				(x1, y11),
				(x1, y12),
				(x2, y21),
				(x2, y22)
			]
			distances = [norm((x, y) .- p) for p in points]
			t = argmin(distances) # time to collision on circle
			xc, yc = points[t]    # collision point on circle

			V = unitdir |> collect
			N = [xc, yc]
			R = V - 2 * ((dot(N, V) / dot(N, N)) * N)
			R /= norm(R)
			return t, (xc, yc), R |> Tuple, true
		end
	elseif δi == 2
		m = (bot - y) / (xx - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with bottom
			return t, (xx, bot), (u, -v), false
		else     # collision with circle
			x1  = (-m*b + √Δ) / (m^2 + 1)
			x2  = (-m*b - √Δ) / (m^2 + 1)
			y11 =  sqrt(r^2 - x1^2)
			y12 = -sqrt(r^2 - x1^2)
			y21 =  sqrt(r^2 - x2^2)
			y22 = -sqrt(r^2 - x2^2)
			points = [
				(x1, y11),
				(x1, y12),
				(x2, y21),
				(x2, y22)
			]
			distances = [norm((x, y) .- p) for p in points]
			t = argmin(distances) # time to collision on circle
			xc, yc = points[t]    # collision point on circle

			V = unitdir |> collect
			N = [xc, yc]
			R = V - 2 * ((dot(N, V) / dot(N, N)) * N)
			R /= norm(R)
			return t, (xc, yc), R |> Tuple, true
		end
	elseif δi == 3
		m = (yy - y) / (left - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with left
			return t, (left, yy), (-u, v), false
		else     # collision with circle
			x1  = (-m*b + √Δ) / (m^2 + 1)
			x2  = (-m*b - √Δ) / (m^2 + 1)
			y11 =  sqrt(r^2 - x1^2)
			y12 = -sqrt(r^2 - x1^2)
			y21 =  sqrt(r^2 - x2^2)
			y22 = -sqrt(r^2 - x2^2)
			points = [
				(x1, y11),
				(x1, y12),
				(x2, y21),
				(x2, y22)
			]
			distances = [norm((x, y) .- p) for p in points]
			t = argmin(distances) # time to collision on circle
			xc, yc = points[t]    # collision point on circle

			V = unitdir |> collect
			N = [xc, yc]
			R = V - 2 * ((dot(N, V) / dot(N, N)) * N)
			R /= norm(R)
			return t, (xc, yc), R |> Tuple, true
		end
	elseif δi == 4
		m = (yy - y) / (right - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with right
			return t, (right, yy), (-u, v), false
		else     # collision with circle
			x1  = (-m*b + √Δ) / (m^2 + 1)
			x2  = (-m*b - √Δ) / (m^2 + 1)
			y11 =  sqrt(r^2 - x1^2)
			y12 = -sqrt(r^2 - x1^2)
			y21 =  sqrt(r^2 - x2^2)
			y22 = -sqrt(r^2 - x2^2)
			points = [
				(x1, y11),
				(x1, y12),
				(x2, y21),
				(x2, y22)
			]
			distances = [norm((x, y) .- p) for p in points]
			t = argmin(distances) # time to collision on circle
			xc, yc = points[t]    # collision point on circle

			V = unitdir |> collect
			N = [xc, yc]
			R = V - 2 * ((dot(N, V) / dot(N, N)) * N)
			R /= norm(R)
			return t, (xc, yc), R |> Tuple, true
		end
	end
end

function billiardtrajectory(
	B::Billiard, 
	δBxy::Tuple{Float64, Float64}, 
	unitdir::Tuple{Float64, Float64},
	steps::Int64
)
	traj = [δBxy]
	uds = [unitdir]
	circle_hits = Int64[]

	for i in 1:steps
		t, newδBxy, newud, circ_hit = nextcollision(B, traj[end], uds[end])

		if circ_hit
			push!(circle_hits, i)
		end

		push!(traj, newδBxy)
		push!(uds, newud)
	end

	return traj, circle_hits
end

function gaps(v::Vector{Int64})
	if length(v) == 0
		return []
	end
	g = Vector{Int64}(undef, length(v)-1)
	for i in 2:length(v)
		g[i - 1] = v[i] - v[i-1]
	end
	return g
end

function randomu0ud(b::Billiard)
	side = rand(1:4)
	if side == 1     # top
		x = (rand() * (last(b.xrange) - first(b.xrange))) - last(b.xrange)
		y = last(b.yrange) |> Float64

		θ = rand() * 2π
		u = cos(θ)
		v = -abs(sin(θ))
	elseif side == 2 # right
		x = last(b.xrange) |> Float64
		y = (rand() * (last(b.yrange) - first(b.yrange))) - last(b.yrange)

		θ = rand() * 2π
		u = -abs(cos(θ))
		v = sin(θ)
	elseif side == 3 # bottom
		x = (rand() * (last(b.xrange) - first(b.xrange))) - last(b.xrange)
		y = first(b.yrange) |> Float64

		θ = rand() * 2π
		u = cos(θ)
		v = abs(sin(θ))
	else             # left
		x = first(b.xrange) |> Float64
		y = (rand() * (last(b.yrange) - first(b.yrange))) - last(b.yrange)

		θ = rand() * 2π
		u = abs(cos(θ))
		v = sin(θ)
	end

	u0 = (x, y)
	ud = (u, v)
	return u0, ud
end


#=
# create circle
xc = 1.5*cos.(0:0.01:2π)
yc = 1.5*sin.(0:0.01:2π)

# create billiard
b = Billiard(-3:3, -3:3, 1.5)

u0 = (-3.0, 1.0)
ud = (1.0, 1.0) ./ norm((1.0, 1.0))
sN = 25
t1 = billiardtrajectory(b, u0, ud, sN)

trajectories = [billiardtrajectory(b, u0 .+ (rand()/100, rand()/100), ud, sN) for _ in 1:100]

anim = @animate for i in 1:sN
	# plot circle
	plot(xc, yc, label=false, xlim=[-3,3], ylim=[-3,3], aspect_ratio=:equal, lc=:black)
	
	points = getindex.(trajectories, i)
	scatter!(points, c=:blue, label=false)
	scatter!(t1[i], c=:red, label=false)

end
gif(anim, "plots/exercise8_4.gif", fps=5)
=#

nBs = []
rs = range(0.5, 2.5; length=100)
samples = 100
sample_length = 1000

for r in rs
	b = Billiard(-3:3, -3:3, r)

	allgaps = []
	for i in 1:samples
		u0, ud = randomu0ud(b)
		tr, ch = billiardtrajectory(b, u0, ud, sample_length)
		push!(allgaps, gaps(ch)...)
	end

	nB = allgaps |> mean
	push!(nBs, nB)
end

plot(
	rs, 
	1 .+ (12 ./ (rs .* π)), 
	ls=:dash,
	lc=:black, 
	label="Theoretical (Kac's Lemma)",
	xlabel="\$ r\$",
	ylabel="Mean return",
	title="Sinai Billiard"
)
scatter!(
	rs, 
	nBs, 
	label="Observed ($sample_length collisions)", 
	mc=:indianred, 
	msw=0,
	ms=2
)
savefig("plots/exercise8_4.png")
