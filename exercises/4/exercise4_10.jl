using Plots; ENV["GKSwstype"] = "100"

function orbit(map, x0, N; p=[2.75])
	xs = [x0]
	for _ in 2:N
		push!(xs, map(xs[end]; p=p))
	end
	return xs
end

function orbit_diagram(map, x0, prange; N=10, t=6, k=1)
	points = []
	for p in prange
		orbitp = orbit(map, x0, N; p=[p])[t:end]
		period = unique(x->round(x,digits=k), orbitp)
		for point in period
			push!(points, (p, point))
		end
	end
	return first.(points), last.(points)
end

function logistic_map(x; p=[2.75])
	r = p[1]
	return r * x * (1 - x)
end

let x0 = 0.4
	scatter(
		orbit_diagram(logistic_map, x0, 2.5:0.0001:4.0; N=1000, t=900, k=2),
		xlim=[2.5,4.0], xtick=2.5:0.5:4.0,
		ylim=[0.0,1.0], ytick=0.0:0.2:1.0,
		msw=0, mc=:black, ms=0.5, mα=0.2,
		label=false,
		title="Orbit Diagram of Logistic Map",
		xlabel="\$r\$", ylabel="\$x\$"
	)
	savefig("plots/exercise4_10.png")
end

let ts = 970:1001
	anim = @animate for t in ts
		println(t)
		scatter(
			orbit_diagram(logistic_map, x0, 2.5:0.0001:4.0; N=1000, t=t, k=2),
			xlim=[2.5,4.0], xtick=2.5:0.5:4.0,
			ylim=[0.0,1.0], ytick=0.0:0.2:1.0,
			msw=0, mc=:black, ms=0.5, mα=0.2,
			label="\$t=$t\$",
			title="Orbit Diagram of Logistic Map",
			xlabel="\$r\$", ylabel="\$x\$"
		)
	end
	gif(anim, "plots/exercise4_10.gif", fps=5)
end
