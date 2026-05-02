using LinearAlgebra
using Plots

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
			return t, (xx, top), (u, -v)
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
			return t, (xc, yc), R |> Tuple
		end
	elseif δi == 2
		m = (bot - y) / (xx - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with bottom
			return t, (xx, bot), (u, -v)
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
			return t, (xc, yc), R |> Tuple
		end
	elseif δi == 3
		m = (yy - y) / (left - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with left
			return t, (left, yy), (-u, v)
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
			return t, (xc, yc), R |> Tuple
		end
	elseif δi == 4
		m = (yy - y) / (right - x)
		b = y - m * x
		Δ = (m^2 * b^2) - ((m^2 + 1) * (b^2 - r^2))

		if x^2 + y^2 ≈ r^2 || Δ < 0 # collision with right
			return t, (right, yy), (-u, v)
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
			return t, (xc, yc), R |> Tuple
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

	for _ in 1:steps
		t, newδBxy, newud = nextcollision(B, traj[end], uds[end])
		push!(traj, newδBxy)
		push!(uds, newud)
	end

	return traj
end


# create circle
xc = 1.5*cos.(0:0.01:2π)
yc = 1.5*sin.(0:0.01:2π)

# create billiard
b = Billiard(-3:3, -3:3, 1.5)

u0 = (-3.0, 1.0)
ud = (1.0, 1.0) ./ norm((1.0, 1.0))
sN = 25
t1 = billiardtrajectory(b, u0, ud, sN)

#plot!(t1, label=false)

trajectories = [billiardtrajectory(b, u0 .+ (rand()/100, rand()/100), ud, sN) for _ in 1:100]

anim = @animate for i in 1:sN
	# plot circle
	plot(xc, yc, label=false, xlim=[-3,3], ylim=[-3,3], aspect_ratio=:equal, lc=:black)
	
	points = getindex.(trajectories, i)
	scatter!(points, c=:blue, label=false)
	scatter!(t1[i], c=:red, label=false)

end
gif(anim, "plots/exercise8_4.gif", fps=1)
