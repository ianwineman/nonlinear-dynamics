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

	ttop = (top - y) / v   # collision time to y = top
	tbot = (bot - y) / v   # collision time to y = bot
	tlef = (left - x) / u  # collision time to x = left
	trig = (right - x) / u # collision time to x = right

	coltimes = [ttop, tbot, tlef, trig]
	t = minimum(filter(>(0), coltimes))
	δi = findfirst(==(t), coltimes)

	xx = x + (t * u)
	yy = y + (t * v)
	if δi == 1     # collision with top
		return t, (xx, top), (u, -v)
	elseif δi == 2 # collision with bottom
		return t, (xx, bot), (u, -v)
	elseif δi == 3 # collision with left
		return t, (left, yy), (-u, v)
	elseif δi == 4 # collision with right
		return t, (right, yy), (-u, v)
	end
end

b = Billiard(-3:3, -3:3, 1.5)

traj = [(-3.0, 1.0)]
ud = (1.0, 1.0) ./ norm((1.0, 1.0))

for _ in 1:10
	_, newxy, newud = nextcollision(b, traj[end], ud)
	global ud = newud
	push!(traj, newxy)
end

plot(traj, xlim=[-3,3], ylim=[-3,3], aspect_ratio=:equal, label="b1")

dx, dy = rand()/10, rand()/10
traj = [(-3.0 + dx, 1.0 + dy)]
ud = (1.0, 1.0) ./ norm((1.0, 1.0))

for _ in 1:10
	_, newxy, newud = nextcollision(b, traj[end], ud)
	global ud = newud
	push!(traj, newxy)
end

plot!(traj, label="($dx, $dy)")
