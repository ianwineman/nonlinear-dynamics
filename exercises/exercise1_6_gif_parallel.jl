using Plots

function standard_map(θ, v; k)
	return (θ + v + k*sin(θ), v + k*sin(θ)) .% 2π
end

function iterate_until_orbit(θ, v; k)
	last = (θ, v)
	for step=1:1000
		s = standard_map(last...,k=k)
		last = s
		if isapprox(s[1],θ,atol=0.1) && isapprox(s[2],v,atol=0.1)
			return step
		else
			if step == 1000
				return 1000
			end
		end
	end
end

function orbits_of_points(points; k)
	orbits = zeros(Int,length(points))
	Threads.@threads for i=1:length(points)
		orbits[i] = iterate_until_orbit(points[i]...,k=k)
	end
	return orbits
end


anim = @animate for k=0:0.01:5
	θs = 0.0:0.1:2π
	vs = 0.0:0.1:2π

	points = [(θ,v) for θ=θs for v=vs]
	orbits = orbits_of_points(points,k=k)

	heatmap(
		θs, vs, reshape(orbits, length(θs), length(vs))', 
		c=cgrad(:default, [0.01, 0.02, 0.03, 0.1, 0.2, 0.3, 1.0]),
		title="Orbits of Standard Map for k=$k",
		xguide="v",
		yguide="θ"
	)

    println("k=$k")
end

gif(anim, "plots/exercise1_6_parallel.gif", fps=5)

# export GKSwstype=nul
# https://github.com/JuliaPlots/Plots.jl/issues/3664#issuecomment-887365869


## BENCHMARKING
## using θs = 0.0:0.1:2π; vs = 0.0:0.1:2π; k=0:0.01:5
#      with multi-threaded orbits_of_points function
#      tested on Macbook
#      n-threads | time
#      1         | 53.97s user 1.92s system 107% cpu 52.199 total
#      2         | 55.95s user 1.67s system 153% cpu 37.643 total
#      4         | 59.04s user 1.67s system 198% cpu 30.645 total
#      8         | 65.30s user 1.92s system 257% cpu 26.111 total
#
## using θs = 0.0:0.01:2π; vs = 0.0:0.01:2π; k=0:0.01:5
#      with multi-threaded orbits_of_points function
#      tested on c7g.16xlarge EC2 instance (64 cores)
#      64 threads: 3m22.822s