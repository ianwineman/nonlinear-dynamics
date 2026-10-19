using LinearAlgebra, Plots, Statistics, StatsBase

function ER(n, p)
	A = Int.(Matrix(LowerTriangular(rand(n, n) .<= p)))
	A .+= transpose(A)
	A .-= Diagonal(A)
end

@enum State susceptible infectious recovered

mutable struct Agent
	state::State
	days_infected::Int
	location::Int
end

function init(n, p, N)
	A = ER(n, p)
	node_degrees = [sum(A[i, :]) for i in 1:n]
	agents = reduce(vcat, [[Agent(susceptible, 0, i) for _ in 1:N] for i in 1:n])
	agents[rand(1:(n * N))].state = infectious

	return A, node_degrees, agents
end

"modifies `agents`"
function step!(A, node_degrees, agents, p, τ, m)
	# recover
	for a in agents
		if a.days_infected > p
			a.state = recovered
			a.days_infected = 0
		end
	end

	# infect neighbors
	for a in agents
		if a.state == infectious
			susceptible_neighbors = [b for b in agents if (b.location == a.location) && (b.state == susceptible)]
			for b in susceptible_neighbors
				if τ > rand()
					b.state = infectious
					b.days_infected = 0
				end
			end
		end
	end

	# travel
	for a in agents
		if m > rand()
			# agent travels
			potential_destinations = findall(==(1), A[a.location, :])
			destination_weights = node_degrees[potential_destinations]
			destination = sample(potential_destinations, Weights(destination_weights))
			a.location = destination
		end
	end

	# update days infected
	for a in agents
		a.days_infected += 1
	end
end

function infectiouscount(agents)
	count(a -> a.state == infectious, agents)
end

trials = 100
steps = 50
d = 0.2
p = 14
τ = 0.3
m = 0.1
Is = []
for i in 1:trials
	i |> println

	A, node_degrees, agents = init(100, d, 10)

	I = zeros(Int, steps)
	I[1] = 1
	for j in 2:steps
		step!(A, node_degrees, agents, p, τ, m)
		I[j] = infectiouscount(agents)
	end
	push!(Is, I)
end
Is = Is |> stack
means = [mean(Is[i, :]) for i in 1:size(Is)[1]]
σs = [std(Is[i, :]) for i in 1:size(Is)[1]]

plot(
	1:steps, 
	means, 
	ribbon=σs, 
	fillalpha=0.3, 
	lw=2,
	xlim=[0, steps],
	label="\$ τ=$τ, m=$m\$",
	xlabel="\$ t\$",
)

trials = 100
steps = 50
d = 0.2
p = 14
τ = 0.15
m = 0.05
Is = []
for i in 1:trials
	i |> println

	A, node_degrees, agents = init(100, d, 10)
	
	I = zeros(Int, steps)
	I[1] = 1
	for j in 2:steps
		step!(A, node_degrees, agents, p, τ, m)
		I[j] = infectiouscount(agents)
	end
	push!(Is, I)
end
Is = Is |> stack
means = [mean(Is[i, :]) for i in 1:size(Is)[1]]
σs = [std(Is[i, :]) for i in 1:size(Is)[1]]

plot!(
	1:steps, 
	means, 
	ribbon=σs, 
	fillalpha=0.3, 
	lw=2,
	xlim=[0, steps],
	label="\$ τ=$τ, m=$m\$",
)
savefig("plots/10.13.png")
