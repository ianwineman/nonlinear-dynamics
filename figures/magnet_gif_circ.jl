#=
# Generate gif of moving magnets path
using Plots; ENV["GKSwstype"] = "100"

θ = range(0, 2π; length=100)
anim = @animate for i in 1:100
	plot(cos.(θ), sin.(θ), aspect_ratio=:equal, lc=:black, lw=1.5, label=false, xlim=[-5, 5], ylim=[-5, 5], xtick=[-5, 0, 5], ytick=[-5, 0, 5])
	scatter!([-0.5, -0.5], [sqrt(3)/2, -sqrt(3)/2], mc=:blue, ms=6, label=false)
	plot!(sqrt(2) .* cos.(θ) .- (sqrt(2) - 1) , sqrt(2) .* sin.(θ), lc=:black, lw=1.5, ls=:dash, label=false)
	scatter!([sqrt(2) * cos(θ[i]) - (sqrt(2) - 1)], [sqrt(2) * sin(θ[i])], mc=:red, ms=6, label=false)
end
gif(anim, "plots/mp/mp_circ.gif", fps=30)
=#

#=
# show sample trajectory
using DifferentialEquations
using Plots, Measures; ENV["GKSwstype"] = "100"

function magnetic_pendulum!(du, u, p , t)
	ω, q, γ, d, x1, y1 = p
	x, y, vx, vy = u

	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2

	du[1] = vx
	du[2] = vy
	du[3] = -(ω^2)*x - q*vx - ((γ*(x-x1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(x-x2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(x-x3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
	du[4] = -(ω^2)*y - q*vy - ((γ*(y-y1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(y-y2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(y-y3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
end

function trajectory(system, u, tspan, p)
	prob = ODEProblem(system, u, tspan, p)
	solu = solve(prob)
	return solu
end

anim = @animate for tend in 0.0:0.25:100.0
	println(tend)
	tr = trajectory(magnetic_pendulum!, [3, -3, 0.0, 0.0], (0.0, tend), [1.0, 0.2, 1.0, 0.3, 1, 0])
	p1 = plot(
		tr, 
		idxs=(1,2), 
		label="tr", 
		xlabel="x", 
		ylabel="y", 
		xlim=[-5, 5], 
		ylim=[-5, 5], 
		xtick=[-5, 0, 1, 5], 
		ytick=[-5, 0, 5],
		aspect_ratio=:equal
	)
	scatter!([3], [-3], label="u0")
	scatter!([1], [0], label="m1")
	scatter!([-0.5], [sqrt(3)/2], label="m2")
	scatter!([-0.5], [-sqrt(3)/2], label="m3")
	p2 = plot(
		tr.t, 
		[getindex.(tr.u, 3) getindex.(tr.u, 4)], 
		label=["ẋ" "ẏ"], 
		xlabel="t", 
		ylabel="v",
		xlim=[0, 100]
	)
	plot(p1, p2, size=(800, 400), margin=5mm)
end
gif(anim, "plots/mp/mp_circ_tr.gif", fps=30)
=#

#=
# show basin of attraction for a given m1
using Images
using DifferentialEquations
using Base.Threads
using ColorSchemes, Colors
using ProgressMeter

function magnetic_pendulum!(du, u, p , t)
	ω, q, γ, d, x1, y1 = p
	x, y, vx, vy = u

	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2

	du[1] = vx
	du[2] = vy
	du[3] = -(ω^2)*x - q*vx - ((γ*(x-x1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(x-x2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(x-x3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
	du[4] = -(ω^2)*y - q*vy - ((γ*(y-y1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(y-y2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(y-y3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
end

function trajectory(system, u, tspan, p)
	prob = ODEProblem(system, u, tspan, p)
	solu = solve(prob)
	return solu
end

function closest_magnet(x, y; m1=(1.0, 0.0))
	x1, y1 = m1
	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2
	δ = [
		sqrt((x-x1)^2 + (y-y1)^2), 
		sqrt((x-x2)^2 + (y-y2)^2), 
		sqrt((x-x3)^2 + (y-y3)^2)
	]
	return argmin(δ)
end

function generate_image(size, theme; m1=(1.0, 0.0), warn=false)
	tspan = (0.0, 100.0)
	p = [1.0, 0.2, 1.0, 0.3]
	img = zeros(RGB{Float64}, size...)
	xs, ys = range(-5.0, 5.0; length=size[1]), range(5.0, -5.0; length=size[2])
	colors = [get(colorschemes[theme], i) for i in 0.0:0.5:1.0]

	@showprogress @threads for i in 1:size[1]
		for j in 1:size[2]
			tr = trajectory(magnetic_pendulum!, [xs[i], ys[j], 0.0, 0.0], tspan, vcat(p, m1...))
			if warn && sum(abs.(tr.u[end][[3, 4]])) > 1e-3
				@warn "tolerance = $(sum(abs.(tr.u[end][[3, 4]]))) > 1e-3"
			end
			cm = closest_magnet(tr.u[end][1:2]...; m1=m1)
			img[j, i] = colors[cm]
		end
	end
	return img
end

img = generate_image((100, 100), :PuRd_3)
save("plots/mp/mp_small.png", img)
=#

#=
colors = [
	:bukavu,
	:buda,
	:davos,
	:fes,
	:grayC,
	:hawaii,
	:imola,
	:lajolla,
	:navia,
	:Blues_3,
	:BrBG_3,
	:BuPu_3,
	:OrRd_5,
	:PuOr_3,
	:PuRd_3
]

for color in colors
	img = generate_image((250, 250), color)
	save("plots/mp/colors/mp_$(string(color)).png", img)
end
=#

# show basins of attraction for m1 as it moves
using Images
using DifferentialEquations
using Base.Threads
using ColorSchemes, Colors
using ProgressMeter

function magnetic_pendulum!(du, u, p , t)
	ω, q, γ, d, x1, y1 = p
	x, y, vx, vy = u

	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2

	du[1] = vx
	du[2] = vy
	du[3] = -(ω^2)*x - q*vx - ((γ*(x-x1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(x-x2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(x-x3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
	du[4] = -(ω^2)*y - q*vy - ((γ*(y-y1))/sqrt((x-x1)^2 + (y-y1)^2 + d^2)^3 + (γ*(y-y2))/sqrt((x-x2)^2 + (y-y2)^2 + d^2)^3 + (γ*(y-y3))/sqrt((x-x3)^2 + (y-y3)^2 + d^2)^3)
end

function trajectory(system, u, tspan, p)
	prob = ODEProblem(system, u, tspan, p)
	solu = solve(prob; save_everystep=false, save_start=false)
	return solu
end

function closest_magnet(x, y; m1=(1.0, 0.0))
	x1, y1 = m1
	x2, y2 = -0.5, sqrt(3)/2
	x3, y3 = -0.5, -sqrt(3)/2
	δ = [
		sqrt((x-x1)^2 + (y-y1)^2), 
		sqrt((x-x2)^2 + (y-y2)^2), 
		sqrt((x-x3)^2 + (y-y3)^2)
	]
	return argmin(δ)
end

function generate_image(size, theme; m1=(1.0, 0.0), warn=false)
	tspan = (0.0, 100.0)
	p = [1.0, 0.2, 1.0, 0.3]
	img = zeros(RGB{Float64}, size...)
	xs, ys = range(-5.0, 5.0; length=size[1]), range(5.0, -5.0; length=size[2])
	colors = [get(colorschemes[theme], i) for i in 0.0:0.5:1.0]
	pm1 = vcat(p, m1...)

	@threads for i in 1:size[1]
		for j in 1:size[2]
			tr = trajectory(magnetic_pendulum!, [xs[i], ys[j], 0.0, 0.0], tspan, pm1)
			cm = closest_magnet(tr.u[end][1:2]...; m1=m1)
			img[j, i] = colors[cm]
		end
	end
	return img
end

function generate_gif(size, theme; warn=false)
	θs = range(0, 2π; length=100)
	imgs = []
	@showprogress for θ in θs
		m1 = sqrt(2) * cos(θ) - (sqrt(2) - 1), sqrt(2) * sin(θ)
		img = generate_image(size, theme; m1=m1, warn=warn)
		push!(imgs, img)
	end
	save("plots/mp/mp_circ_basins_$(size[1])x$(size[2])_$(string(theme)).gif", cat(imgs..., dims=3); fps=10)
	println("plots/mp/mp_circ_basins_$(size[1])x$(size[2])_$(string(theme)).gif")
end

generate_gif((1_000, 1_000), :BrBG_3)
