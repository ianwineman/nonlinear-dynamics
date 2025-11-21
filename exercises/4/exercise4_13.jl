using DifferentialEquations
using Plots
#using Base.Threads

function gissinger!(du, u, p, t)
	μ, ν, γ = p
	x, y, z = u

	du[1] = μ*x - y*z
	du[2] = -ν*y + x*z
	du[3] = γ - z + x*y
end

function neg_to_pos(vec)
	indices = []
	neg = findall(x->x<0, vec)
	pos = findall(x->x>0, vec)

	for i in neg
		if (i + 1) in pos
			push!(indices, (i, i+1))
		end
	end
	return indices
end

function psos(system, u0, p, tspan; plane=1, transient_ratio=0.75)
	prob = ODEProblem(system, u0, tspan, p)
	solu = solve(prob, dtmax=0.01)

	cross_plane = neg_to_pos(getindex.(solu.u, plane))
	ts = [solu.t[pair[1]] + (solu.t[pair[2]]-solu.t[pair[1]])/2 for pair in cross_plane]
	return [solu(t) for t in filter(x->x>tspan[2]*transient_ratio, ts)]
end

function orbit_points(p_index, p_range, p_default, u0)
	orbit_points = []
	for param in p_range # how to parallelize ???
		println("$param ∈ $p_range")
		p_default[p_index] = param
		psos_points = psos(gissinger!, u0, p_default, (0.0,1000.0); plane=3, transient_ratio=0.9)
		for p in psos_points
			push!(orbit_points, (param, p[1]))
		end
	end
	return vcat(orbit_points...)
end

op = orbit_points(1, 0.1:0.00001:0.15, [0.119, 0.1, 0.9], [3.0, 3.0, 3.0])

scatter(
   first.(op), last.(op),
   label=false, xlabel="\$μ\$", ylabel="\$x\$",
   title="Orbit diagram of Gissinger system",
   msw=0, mc=:black, ms=0.5, mα=0.5
)
savefig("plots/exercise4_13_mu.png")

op = orbit_points(2, 0.08:0.00001:0.18, [0.119, 0.1, 0.9], [3.0, 3.0, 3.0])

scatter(
   first.(op), last.(op),
   label=false, xlabel="\$ν\$", ylabel="\$x\$",
   title="Orbit diagram of Gissinger system",
   msw=0, mc=:black, ms=0.5, mα=0.5
)
savefig("plots/exercise4_13_nu.png")
