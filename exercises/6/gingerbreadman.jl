using Plots, ProgressMeter

# https://en.wikipedia.org/wiki/Gingerbreadman_map
function gingerbreadmanmap(u)
	x, y = u
	return [1 - y + abs(x), x]
end

xs = -10.0:0.01:10
ys = -10.0:0.01:10

points = [[x,y] for x=xs for y=ys]
orbits = zeros(Int, length(points))

@showprogress Threads.@threads for (i, p) in collect(enumerate(points))
	steps = [p]
	for step=1:1000
		s = gingerbreadmanmap(steps[end])
		push!(steps,s)
		if isapprox(s[1],p[1],atol=0.1) && isapprox(s[2],p[2],atol=0.1)
			orbits[i] = step
			break
		else
			if step == 1000
				orbits[i] = 1000
			end
		end
	end
end

heatmap(
   xs, ys, reshape(orbits, length(xs), length(ys))',
   c=cgrad(:devon),
   title="Gingerbread man map",
   xguide="\$x\$",
   yguide="\$y\$",
   colorbar_title="Orbit length"
)
savefig("plots/gingerbreadman.png")
