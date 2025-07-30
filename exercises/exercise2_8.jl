using Plots, Measures, DifferentialEquations, LinearAlgebra

function brusselator!(du, u, p, t)
   a, b = p
   x, y = u
   du[1] = 1 - (b+1)*x + a*(x^2)*y
   du[2] = b*x - a*(x^2)*y
end

function angle(u,w,u_fp,w_fp)
	if u >= u_fp && w >= w_fp # quadrant I
		return atan(abs(w-w_fp)/abs(u-u_fp))
	elseif u <= u_fp && w >= w_fp # quadrant II
		return π - atan(abs(w-w_fp)/abs(u-u_fp))
	elseif u <= u_fp && w <= w_fp # quadrant III
		return π + atan(abs(w-w_fp)/abs(u-u_fp))
	elseif u >= u_fp && w <= w_fp # quadrant IV
		return 2*π - atan(abs(w-w_fp)/abs(u-u_fp))
	end
end

function arclength(sol, ϵ, s1, s2)
	arc = collect(
				Iterators.takewhile(
					x -> (norm(sol.u[1] .- x) > ϵ),
					sol.u[2:end]
				)
			)
	
	arclength = []
	for (i, _) in enumerate(arc)
		if i == length(arc)
			push!(arclength, abs(norm(arc[1] .- arc[i])))
		else
			push!(arclength, abs(norm(arc[i+1] .- arc[i])))
		end
	end
	
	arc_t = sol.t[2:length(arc)+1]
	return (
		[(arc_t[i],sum(arclength[1:i])/sum(arclength[1:end])) for (i,_) in enumerate(arc)], 
		sum(arclength[1:s1])/sum(arclength[1:end]), 
		sum(arclength[1:s2])/sum(arclength[1:end])
	)
end

tmax = 12.9
prob = ODEProblem(brusselator!, [0.5909924611496271,4.843385159662832], (0.0,tmax), [0.3, 1.5])
sol = solve(prob, dtmax=0.01)

sample1, sample2 = 100, 750

p1 = plot(
	sol, idxs=(1,2), 
	label=false,
	xlabel="u",
	ylabel="w",
	title="Brusselator",
	titlefontsize=12,
	left_margin=5mm,
	linecolor=:blue
)
scatter!(
	[1.0], [5.0],
	label=false,
	mc=:white
)
scatter!(
	[1.0], [4.0],
	label=false,
	mc=:white
)
scatter!(
	[sol.u[sample1][1]], [sol.u[sample1][2]],
	label=false,
	mc=:orange
)
scatter!(
	[sol.u[sample2][1]], [sol.u[sample2][2]],
	label=false,
	mc=:red
)

p2 = plot(
	sol.t,
	[angle(u,w,1,5) for (u,w) in sol.u],
	label="\$θ(t)\$",
	xlim=[0,tmax],
	ylim=[0.0,2*π],
	title="Angle protophase about \$ (1,5)\$",
	titlefontsize=12,
	linecolor=:blue
)
scatter!(
	[sol.t[sample1][1]], [angle(sol.u[sample1][1],sol.u[sample1][2],1,5)],
	label=false,
	mc=:orange
)
scatter!(
	[sol.t[sample2][1]], [angle(sol.u[sample2][1],sol.u[sample2][2],1,5)],
	label=false,
	mc=:red
)

p3 = plot(
	sol.t,
	[angle(u,w,1,4) for (u,w) in sol.u],
	label="\$θ(t)\$",
	xlim=[0,tmax],
	ylim=[0.0,2*π],
	title="Angle protophase about \$ (1,4)\$",
	titlefontsize=12,
	linecolor=:blue
)
scatter!(
	[sol.t[sample1][1]], [angle(sol.u[sample1][1],sol.u[sample1][2],1,4)],
	label=false,
	mc=:orange
)
scatter!(
	[sol.t[sample2][1]], [angle(sol.u[sample2][1],sol.u[sample2][2],1,4)],
	label=false,
	mc=:red
)

p4 = plot(
	arclength(sol,0.001, sample1, sample2)[1],
	title="Arclength protophase",
	titlefontsize=12,
	label="\$ l(t)\$",
	ytick=[0,0.5,1],
	linecolor=:blue
)
scatter!(
	[sol.t[sample1][1]], [arclength(sol,0.001, sample1, sample2)[2]],
	label=false,
	mc=:orange
)
scatter!(
	[sol.t[sample2][1]], [arclength(sol,0.001, sample1, sample2)[3]],
	label=false,
	mc=:red
)
hline!(
	[0.5], 
	linecolor=:black, 
	linestyle=:dash,
	linealpha=0.25
)

plot(
	p1,p2,p3,p4,
	size=(1200,600),
	top_margin=5mm,
	bottom_margin=5mm,
	layout=@layout [ a{0.25w} [grid(1,2)
							b ]]
)
savefig("plots/exercise2_8.png")
