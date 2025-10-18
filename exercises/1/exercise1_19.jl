using DifferentialEquations
using Plots; ENV["GKSwstype"] = "100"

function lorenz84!(du, u, p, t)
	a, b, F, G = p
	x, y, z = u

	du[1] = a*F - a*x - y^2 - z^2 
	du[2] = x*y - y - b*x*z + G
	du[3] = b*x*y + x*z - z
end

function trajectory(u0,tlim)
	tspan = (0.0,tlim)
	p = [0.25, 4.0, 6.886, 1.337]

	prob = ODEProblem(lorenz84!, u0, tspan, p)
	sol = solve(prob, dtmax=0.1)

	return sol
end

anim = @animate for t=1.0:0.1:25.0
	sol_fixed_point_1 = trajectory([0.0,1.337,0.0], t)
	sol_fixed_point_2 = trajectory([0.0,1.337,0.1], t)
	sol_limit_cycle_1 = trajectory([1.337,1.0,0.0], t)
	sol_limit_cycle_2 = trajectory([1.337,1.0,0.1], t)
	sol_chaotic_att_1 = trajectory([-1.0,0.0,0.0], t)
	sol_chaotic_att_2 = trajectory([-1.0,0.1,0.0], t)

	p1 = plot(
		xlabel="x", ylabel="y", zlabel="z", 
		title="Attracting Fixed Point", 
		legendposition=:topright
	)
	plot!(
		sol_fixed_point_1, idxs=(1,2,3), 
		label="\$u_0=$(sol_fixed_point_1.prob.u0)\$", 
		color=:green
	)
	plot!(
		sol_fixed_point_2, idxs=(1,2,3), 
		label="\$u_0=$(sol_fixed_point_2.prob.u0)\$", 
		color=:purple
	)

	p2 = plot(
		xlabel="x", ylabel="y", zlabel="z", 
		title="Attracting Limit Cycle", 
		legendposition=:topright
	)
	plot!(
		sol_limit_cycle_1, idxs=(1,2,3), 
		label="\$u_0=$(sol_limit_cycle_1.prob.u0)\$", 
		color=:green
	)
	plot!(
		sol_limit_cycle_2, idxs=(1,2,3), 
		label="\$u_0=$(sol_limit_cycle_2.prob.u0)\$", 
		color=:purple
	)

	p3 = plot(
		xlabel="x", ylabel="y", zlabel="z", 
		title="Chaotic Attractor",
		legendposition=:topright
	)
	plot!(
		sol_chaotic_att_1, idxs=(1,2,3), 
		label="\$u_0=$(sol_chaotic_att_1.prob.u0)\$", 
		color=:green
	)
	plot!(
		sol_chaotic_att_2, idxs=(1,2,3), 
		label="\$u_0=$(sol_chaotic_att_2.prob.u0)\$", 
		color=:purple
	)

	plot(p1,p2,p3, layout=(1,3), size=(900,300))
	println("t=$t;")
end

gif(anim, "plots/exercise1_19.gif", fps=20);