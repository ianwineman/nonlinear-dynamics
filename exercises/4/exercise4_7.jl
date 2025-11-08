using Plots, LinearAlgebra

function henon_map(u; p=[1.4, 0.3])
	x, y = u
	a, b = p
	return [1 - a*x^2 + y, b*x]
end

function inbasin(map, u0, attractor; maptol=0.001, atttol=0.1)
	u = u0
	while norm(map(u)) < 5.0 && norm(map(u) .- u) > maptol
		u = map(u)
	end
	if norm(u .- attractor) < atttol
		return true
	else
		return false
	end
end

function getbasin(map, x0s, y0s, attractor)
	u0s   = [[x,y] for x in x0s for y in y0s]
	basin = zeros(Int, length(u0s))

	Threads.@threads for i in 1:length(u0s)
		if inbasin(map, u0s[i], attractor)
			basin[i] = i
		end
	end
	return u0s[filter(x->x!=0,basin)]
end

chaotic_attractor = [(-7+10*sqrt(6.09))/28, (-21+30*sqrt(6.09))/280]
x0s = -1.5:0.01:1.5
y0s = -0.5:0.01:0.5
basin = getbasin(henon_map, x0s, y0s, chaotic_attractor)

scatter(
	first.(basin), last.(basin),
	xlim=[-1.5,1.5], ylim=[-0.5,0.5],
	mc=:beige, msw=0,
	label=false, markershape=:rect,
	backgroundinside=:teal,
	grid=false,
	xlabel="\$x\$", ylabel="\$y\$",
	xtick=[-1.5,1.5], ytick=[-0.5,0.5],
	title="Hénon map"
)
scatter!((chaotic_attractor[1], chaotic_attractor[2]), label="Chaotic attractor", mc=:black)
savefig("plots/exercise4_7.png")
