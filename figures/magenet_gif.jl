using Images
using DifferentialEquations
using Base.Threads
using ColorSchemes, Colors
using ProgressMeter

function magnetic_pendulum!(du, u, p , t)
	ω, q, γ, d = p
	x, y, vx, vy = u

	x1, y1 = 1.0, 0.0
	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2

	du[1] = vx # = ẋ
	du[2] = vy # = ẏ
	du[3] = -(ω^2)*x - q*vx - ((γ*(x-x1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(x-x2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(x-x3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3) # = v̇x
	du[4] = -(ω^2)*y - q*vy - ((γ*(y-y1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(y-y2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(y-y3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3) # = v̇y
end

function trajectory(system, u, tspan, p)
	prob = ODEProblem(system, u, tspan, p)
	solu = solve(prob)
	return solu.u
end

function closest_magnet(x,y)
	x1, y1 = 1.0, 0.0
	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2
	δ = [
		sqrt((x-x1)^2 + (y-y1)^2), 
		sqrt((x-x2)^2 + (y-y2)^2), 
		sqrt((x-x3)^2 + (y-y3)^2)
	]

	return argmin(δ)
end

using Images
using DifferentialEquations
using Base.Threads
using ColorSchemes, Colors

function magnetic_pendulum!(du, u, p , t)
	ω, q, γ, d, x2, y2 = p
	x, y, vx, vy = u

	x1, y1 = 1.0, 0.0
	#x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2

	du[1] = vx # = ẋ
	du[2] = vy # = ẏ
	du[3] = -(ω^2)*x - q*vx - ((γ*(x-x1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(x-x2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(x-x3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3) # = v̇x
	du[4] = -(ω^2)*y - q*vy - ((γ*(y-y1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(y-y2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(y-y3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3) # = v̇y
end

function trajectory(system, u, tspan, p)
	prob = ODEProblem(system, u, tspan, p)
	solu = solve(prob)
	return solu.u
end

function closest_magnet(x, y; m2=(-0.5, sqrt(3)/2))
	x1, y1 = 1.0, 0.0
	x2, y2 = m2
	x3, y3 = -0.5, -sqrt(3)/2
	δ = [
		sqrt((x-x1)^2 + (y-y1)^2), 
		sqrt((x-x2)^2 + (y-y2)^2), 
		sqrt((x-x3)^2 + (y-y3)^2)
	]
	return argmin(δ)
end

function generate_image(size, theme; m2=(-0.5, sqrt(3)/2))
	tspan = (0.0,100.0)
	p = [1.0, 0.2, 1.0, 0.3]
	img = zeros(RGB{Float64}, size...)
	xs, ys = range(-5.0,5.0; length=size[1]), range(-5.0,5.0; length=size[2])
	colors = [get(colorschemes[theme], i) for i in 0.0:0.5:1.0]

	for i in 1:size[1]
		@threads for j in 1:size[2]
			t = trajectory(magnetic_pendulum!, [xs[i], ys[j], 0.0, 0.0], tspan, vcat(p, m2...))
			cm = closest_magnet(t[end][1:2]...; m2=m2)
			img[i, j] = colors[cm]
		end
	end
	return img
end

frames = 100
imgs = []
@showprogress dt=0.5 for (i, m2) in enumerate(zip(range(-0.5, 0.5, frames), range(sqrt(3)/2, -sqrt(3)/2, frames)))
	#println("$i/$frames")
	img = generate_image((400, 400), :lapaz; m2=m2)
	push!(imgs, img)
end
save("plots/mp/mp.gif", cat(imgs..., dims=3); fps=10)
