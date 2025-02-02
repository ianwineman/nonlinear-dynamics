using Plots

function standard_map(θ, v; k)
	return (θ + v + k*sin(θ), v + k*sin(θ)) .% 2π
end

θs = 0.0:0.1:2π
vs = 0.0:0.1:2π

points = [(θ,v) for θ=θs for v=vs]
orbits = []

for p=points
	steps = [p]
	for step=1:1000
		s = standard_map(steps[end]...,k=1)
		push!(steps,s)
		if isapprox(s[1],p[1],atol=0.1) && isapprox(s[2],p[2],atol=0.1)
			push!(orbits,step)
			break
		else
			if step == 1000
				push!(orbits,1000)
			end
		end
	end
end

heatmap(
	θs, vs, reshape(orbits, length(θs), length(vs))', 
	c=cgrad(:default, [0.01, 0.02, 0.03, 0.1, 0.2, 0.3, 1.0]),
	title="Orbits of Standard Map for k=1",
	xguide="v",
	yguide="θ"
)

savefig("plots/exercise1_6.png")