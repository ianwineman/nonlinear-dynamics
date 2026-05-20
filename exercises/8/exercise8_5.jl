using Plots

function logisticmap(x)
	return 4 * x * (1 - x)
end

trajectories = [[rand()] for _ in 1:1_000]
iterations = 10_000

for trajectory in trajectories
	for _ in 1:iterations
		push!(trajectory, logisticmap(trajectory[end]))
	end
end

anim = @animate for n in 0:iterations
	histogram(
		getindex.(trajectories, n + 1),
		label="n = $n", 
		legendposition=:top, 
		bins=range(0, 1; length=25), 
		xlim=[0, 1],
		#ylim=[0, 200],
		normalize=:probability
	)
	#x = range(0, 1; length=1_000)
	#plot!(
	#	x,
	#	1 ./ (π .* sqrt.(x .* (1 .- x)))
	#)
end
gif(anim, "plots/exercise8_5.gif", fps=30)
