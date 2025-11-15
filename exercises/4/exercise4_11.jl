using Plots, Base.Threads

function orbit(map, x0, N; p=[2.75])
	xs = zeros(N)
	xs[1] = x0
	for i in 2:N
		xs[i] = map(xs[i-1]; p=p)
	end
	return xs
end

function orbit_periods(map, x0, prange; N=10, t=6, k=1)
	ps      = collect(prange)
	rs      = zeros(length(prange))
	periods = zeros(Int, length(prange))
	@threads for i in 1:length(ps)
		p = ps[i]
		orbitp = orbit(map, x0, N; p=[p])[t:end]
		period = unique(x->round(x,digits=k), orbitp)
		rs[i] = p
		periods[i] = length(period)
	end
	return collect(zip(rs,periods))
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

ops = orbit_periods(logistic_map, BigFloat(0.4), 2.99:0.000001:3.6; N=100_000, t=99_800, k=6)
b   = bifurcations(ops)
println("δ ≈ $(round((b[end-1]-b[end-2])/(b[end]-b[end-1]), digits=4))")
