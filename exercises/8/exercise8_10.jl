using LinearAlgebra
using Plots
using Statistics, StatsBase

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
	ts = []

	for i in 1:steps
		t, newδBxy, newud, circ_hit = nextcollision(B, traj[end], uds[end])

		if circ_hit
			push!(circle_hits, i)
		end

		push!(ts, t)
		push!(traj, newδBxy)
		push!(uds, newud)
	end

	return traj, circle_hits, ts
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

function recurrencematrix(traj, ϵ)
	R = zeros(Int, length(traj), length(traj))
	for i in 1:length(traj)
		for j in 1:length(traj)
			if i != j && norm(traj[i] .- traj[j]) < ϵ
				R[i, j] = 1
			end
		end
	end
	return R
end

function diagonal(A, i)
	[A[size(A)[1] - a, i + a + 1] for a in 0:(size(A)[1] - i - 1)]
end

function diagonallengths(R)
	diagonal_lengths = []
	for i in 1:(size(R)[1] - 1)
		full_diagonal = diagonal(R, i)
		if sum(full_diagonal) >= 2
			current_diagonal_length = 0

			for j in 1:length(full_diagonal)
				if full_diagonal[j] == 1
					current_diagonal_length += 1
				else
					if current_diagonal_length > 0
						push!(diagonal_lengths, current_diagonal_length)
						current_diagonal_length = 0
					end
				end
			end
		end
	end
	return filter(>(1), diagonal_lengths)
end


#=
b = Billiard(-3:3, -3:3, π/2)
u0, ud = randomu0ud(b)
#tr, _ = billiardtrajectory(b, u0, ud, 10_000)

tr = billiardtrajectoryuntilrecurrence(b, u0, ud, 0.1)
=#
#=
plot(
	tr, 
	xlim=[-3, 3],
    ylim=[-3, 3],
	aspect_ratio=:equal, 
	lw=0.1, 
	label=false,
	grid=false,
	xaxis=false,
	yaxis=false,
	size=(1179, 1179),
	background_color=:black,
	lc=:blue
)
savefig("plots/ergodic_billiard.png")
=#


ϵs = [0.001, 0.01, 0.1]
rs = [1.0, 1.5, 2.0]
sample_size = 100
mean_recurrence_times = zeros(length(ϵs), length(rs))
entropy_recurrence_times = zeros(length(ϵs), length(rs))
mean_dlls = zeros(length(ϵs), length(rs))
trajectory_length = 10_000

for (i, ϵ) in enumerate(ϵs)
	for (j, r) in enumerate(rs)
		b = Billiard(-3:3, -3:3, r)
		reccurrence_times = []
		reccurrence_probs = []
		dlls = []

		for s in 1:sample_size
			u0, ud = randomu0ud(b)
			tr, _, _ = billiardtrajectory(b, u0, ud, trajectory_length)

			recurrence_time = findfirst(x -> norm(x .- tr[1]) < ϵ, tr[2:end])
			recurrence_prob = length(findall(x -> norm(x .- tr[1]) < ϵ, tr[2:end])) / trajectory_length

			if !isnothing(recurrence_time)
				push!(reccurrence_times, recurrence_time)
			end
			if recurrence_prob > 0.0
				push!(reccurrence_probs, recurrence_prob)
			end

			push!(dlls, diagonallengths(recurrencematrix(tr, ϵ))...)
		end

		mean_recurrence_times[i, j] = mean(reccurrence_times)
		entropy_recurrence_times[i, j] = -log2(mean(reccurrence_probs)) + 1 # https://www.pik-potsdam.de/members/kurths/publikationen/2009/phys-lett-a-murilo.pdf
		mean_dlls[i, j] = mean(dlls)
		println("(ϵ, r) = ($ϵ, $r) MRT: $(mean(reccurrence_times)) ERT: $(-log2(mean(reccurrence_probs)) + 1) MDLL: $(mean(dlls))")
	end
end

p1 = heatmap(
	rs,
	ϵs,
	mean_recurrence_times,
	xlabel="\$ R\$, Disk radius",
	ylabel="\$ ϵ\$, Recurrence threshold",
	xtick=rs,
	ytick=ϵs,
	c=:devon,
	yscale=:log10,
	cbar_title="\n\$ ⟨r⟩\$, Mean reccurence time \$ (n = $sample_size) \$",
	leftmargin=5Plots.mm,
	rightmargin=5Plots.mm
)

p2 = heatmap(
	rs,
	ϵs,
	entropy_recurrence_times,
	xlabel="\$ R\$, Disk radius",
	ylabel="\$ ϵ\$, Recurrence threshold",
	xtick=rs,
	ytick=ϵs,
	c=:devon,
	yscale=:log10,
	cbar_title="\n\$ H_r\$, Entropy of reccurence times \$ (n = $sample_size) \$",
	leftmargin=5Plots.mm,
	rightmargin=5Plots.mm
)

p3 = heatmap(
	rs,
	ϵs,
	mean_dlls,
	xlabel="\$ R\$, Disk radius",
	ylabel="\$ ϵ\$, Recurrence threshold",
	xtick=rs,
	ytick=ϵs,
	c=:devon,
	yscale=:log10,
	cbar_title="\n\$ ⟨ℓ⟩\$, Mean diagonal line length \$ (n = $sample_size) \$",
	leftmargin=5Plots.mm,
	rightmargin=5Plots.mm
)

plot(p1, p2, p3, layout=(2, 2), size=(1200, 800))
savefig("plots/exercise8_10.png")

#=
b = Billiard(-3:3, -3:3, 2.5)
u0, ud = randomu0ud(b)
traj, _, _ = billiardtrajectory(b, u0, ud, 10000)
R = recurrencematrix(traj, 0.1)
heatmap(R, c=:devon)

diagonallengths(R) |> mean
=#
