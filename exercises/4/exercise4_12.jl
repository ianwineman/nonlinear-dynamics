using Plots

function orbit(map, x0, N; p=[1.0])
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

let x0 = 0.4
	function sin_map(x; p=[1.0])
		r = p[1]
		return r*sin(π*x)
	end

	scatter(
		orbit_diagram(sin_map, x0, 0.0:0.0001:1.0; N=1000, t=900, k=2),
		xlim=[0.0,1.0], xtick=0.0:0.5:1.0,
		ylim=[0.0,1.0], ytick=0.0:0.5:1.0,
		msw=0, mc=:blue, ms=0.5, mα=0.2,
		label=false,
		title="Orbit Diagram of \$x_{n+1}=rsin(πx_n)\$",
		xlabel="\$r\$", ylabel="\$x\$"
	)
	savefig("plots/exercise4_12_sin.png")
end

let x0 = 0.4
	function cos_map(x; p=[1.0])
		r = p[1]
		return r*cos(π*x)
	end

	scatter(
		orbit_diagram(cos_map, x0, 0.0:0.001:3.0; N=1000, t=900, k=2),
		xlim=[0.0,3.0], xtick=0.0:0.5:3.0,
		ylim=[-3.0,3.0], ytick=-3.0:0.5:3.0,
		msw=0, mc=:blue, ms=0.5, mα=0.2,
		label=false,
		title="Orbit Diagram of \$x_{n+1}=rcos(πx_n)\$",
		xlabel="\$r\$", ylabel="\$x\$"
	)
	savefig("plots/exercise4_12_cos.png")
end

let x0 = 0.4
	function xe_map(x; p=[1.0])
		r = p[1]
		return x*exp(-r*(1-x))
	end

	scatter(
		orbit_diagram(xe_map, x0, -3.0:0.001:0.0; N=1000, t=900, k=2),
		xlim=[-3.0, 0.0], xtick=-3.0:0.5:0.0,
		ylim=[0.0,2.5], ytick=0.0:0.5:2.5,
		msw=0, mc=:blue, ms=0.5, mα=0.2,
		label=false,
		title="Orbit Diagram of \$x_{n+1}=x_ne^{-r(1-x_n)}\$",
		xlabel="\$r\$", ylabel="\$x\$"
	)
	savefig("plots/exercise4_12_xe.png")
end

let x0 = 0.4
	function e_map(x; p=[1.0])
		r = p[1]
		return exp(-r*x)
	end

	scatter(
		orbit_diagram(e_map, x0, -1.0:0.001:5.0; N=1000, t=900, k=2),
		xlim=[-1.0, 5.0], xtick=-1.0:0.5:5.0,
		ylim=[0.0,2.5], ytick=0.0:0.5:2.5,
		msw=0, mc=:blue, ms=0.5, mα=0.2,
		label=false,
		title="Orbit Diagram of \$x_{n+1}=e^{-rx_n}\$",
		xlabel="\$r\$", ylabel="\$x\$"
	)
	savefig("plots/exercise4_12_e.png")
end
