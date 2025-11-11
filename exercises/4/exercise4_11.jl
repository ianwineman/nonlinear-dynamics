using Plots

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

function orbit_periods(map, x0, prange; N=10, t=6, k=1)
	periods = []
	for p in prange
		orbitp = orbit(map, x0, N; p=[p])[t:end]
		period = unique(x->round(x,digits=k), orbitp)
		push!(periods, (p,length(period)))
	end
	return periods
end

function bifurcations(orb_periods)
	k = 0
	rks = [] # rks[k] = r_k
	for (r,p) in orb_periods
		if p == 2^(k + 1)
			push!(rks, r)
			k += 1
		end
	end
	rks
end

function logistic_map(x; p=[2.75])
	r = p[1]
	return r * x * (1 - x)
end

#orbit_diagram(logistic_map, 0.4, 2.5:0.0001:3.6; N=10_000, t=1_000, k=3)
#orbit_periods(logistic_map, 0.4, 2.5:0.0001:3.6; N=10_000, t=1_000, k=3)

#sort(unique(last.(orbit_periods(logistic_map, 0.4, 2.5:0.0001:4.0; N=10_000, t=1_000, k=3))))
