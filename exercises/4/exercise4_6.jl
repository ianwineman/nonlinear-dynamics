using DifferentialEquations, Plots

function system!(du, u, p, t)
	p = p[1]
    x, y = u
    du[1] = y
    du[2] = -p*y + x - x^2 + x*y
end

let 
	u0s = [[0.01, 0.0], [0.06, 0.0], [0.275, 0.0], [0.48, 0.0], [0.85, 0.0]]
	ps = [0.86, 0.87, 0.91, 0.95, 1.0]
	ts = [17.0, 15.0, 15.0, 15.0, 15.0]
	plt = scatter(
		(0,0), 
		label="Saddle point", 
		mc=:black,
		xlabel="\$x\$", ylabel="\$y\$",
		title="\$ẋ=y\$ ; \$ẏ=-py + x - x^2 + xy\$",
		aspectratio=:equal
	)
	for (i, p) in enumerate(ps)
		prob = ODEProblem(system!, u0s[i], (0.0,ts[i]), [p])
		sol = solve(prob)
		plot!(sol, idxs=(1,2), label="\$p=$p\$")
	end
	savefig("plots/exercise4_6.png")
end
