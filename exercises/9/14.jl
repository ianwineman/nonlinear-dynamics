using DifferentialEquations
using Plots

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
	x, y, z = u

	du[1] = σ * (y - x)
	du[2] = -x * z + ρ * x - y
	du[3] = x * y - β * z
end

function coupledlorenz63!(du, u, p, t)
	σ, ρ, r, β = p
	x1, x2, x3, y1, y2, y3 = u

	du[1] = σ * (x2 - x1)
	du[2] = -x1 * x3 + ρ * x1 - x2
	du[3] = x1 * x2 - β * x3
	du[4] = σ * (y2 - y1)
	du[5] = -y1 * y3 + r * x1 - y2
	du[6] = y1 * y2 - β * y3
end


rs = 1.0:50.0
Es = zeros(length(rs))
for (i, r) in enumerate(rs)
	"$i/$(length(rs))" |> println
	prob = ODEProblem(coupledlorenz63!, [20.0, 20.0, 20.0, 20.0, 20.0, 20.0], (0.0, 100.0), [10.0, 28.0, r, 8/3])
	solv = solve(prob, dtmax=0.01)
	x1 = getindex.(solv.u, 1)
	y1 = getindex.(solv.u, 4)
	Es[i] = sum((x1 .- y1) .^2) / length(x1)
end

plot(rs, Es, label="\$ E(r)\$", xlabel="\$ r\$", title="\$ r=28\$")
scatter!([argmin(Es)], [Es[argmin(Es)]], label="argmin")
savefig("plots/9.14.png")
