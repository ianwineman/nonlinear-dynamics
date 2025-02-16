using DifferentialEquations
using LinearAlgebra
using Plots

function duffing!(du, u, p, t)
	d = p
	x, y = u

	du[1] = y
	du[2] = x - x^3 - d*y
end

function trajectory(u0)
	tspan = (0.0,100.0)
	p = 0.2

	prob = ODEProblem(duffing!, u0, tspan, p)
	sol = solve(prob)

	return sol.u[end]
end

function attractor(u0, attractors)
	for a=attractors
		if isapprox(norm(trajectory(u0).-a), 0.0; atol=0.001)
			return findfirst(x -> x == a, attractors)
		end
	end
	return 0
end

function zs_of_points(points, attractors)
	zs = zeros(Float64,length(points))
	Threads.@threads for i=1:length(points)
		zs[i] = attractor(points[i],attractors)
	end
	return zs
end

xs = range(-5,5,length=1000)
ys = range(-5,5,length=1000)
points = [[x,y] for x=xs for y=ys]
zs = zs_of_points(points,[[-1.0,0.0],[1.0,0.0]])

J = [0 1; 1 -0.2]
v1 = eigen(J).vectors[:,1]
v2 = eigen(J).vectors[:,2]
λ1 = eigen(J).values[1]
λ2 = eigen(J).values[2]
m1 = v1[1]/v1[2]
m2 = v2[1]/v2[2]

contour(
   xs,ys,zs,
   fill=true,
   color=[:mediumaquamarine,:mediumpurple],
   colorbar=:none,
   xlabel="x", ylabel="y", title="Duffing oscillator \$(d=0.2)\$",
   xlim=[-5,5], ylim=[-5,5],
   xtick=[-5,-4,-3,-2,-1,0,1,2,3,4,5],
   ytick=[-5,-4,-3,-2,-1,0,1,2,3,4,5],
)
scatter!(
	[-1.0,1.0], [0.0,0.0], 
	ms=4, msw=2, mc=[:mediumpurple,:mediumaquamarine], 
	label=:none
)
plot!(
	range(-5,5,length=100), m1.*range(-5,5,length=100),
	lw=4, lc=:black, label="\$M_s\$"
)
plot!(
	range(-5,5,length=100), m2.*range(-5,5,length=100),
	lw=4, lc=:orange, label="\$M_u\$"
)
quiver!(
	[-2, 2],[m1*-2, m1*2], quiver=([1, -1],[m1, m1*-1]),
	lw=4, color=:black
)
quiver!(
	[0, 0],[0, 0], quiver=([-1, 1],[m2*-1, m2]),
	lw=4, color=:orange
)
scatter!(
	[0.0], [0.0], 
	ms=4, msw=2, mc=:black, 
	label=:none
)

savefig("plots/exercise1_16.png")