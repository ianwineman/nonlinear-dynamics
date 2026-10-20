using LinearAlgebra, Plots, Statistics, StatsBase

function ER(n, p)
	A = Int.(Matrix(LowerTriangular(rand(n, n) .<= p)))
	A .+= transpose(A)
	A .-= Diagonal(A)
end

@enum State susceptible infectious recovered
@enum QState quarantined notquarantined

mutable struct Agent
	state::State
	days_infected::Int
	location::Int
end

mutable struct Node
	degree::Int
	quarantine_state::QState
	days_quarantined::Int
end

function init(n, p, N)
	A = ER(n, p)
	node_degrees = [sum(A[i, :]) for i in 1:n]
	agents = reduce(vcat, [[Agent(susceptible, 0, i) for _ in 1:N] for i in 1:n])
	agents[rand(1:(n * N))].state = infectious

	return A, Node.(node_degrees, [notquarantined for _ in 1:n], zeros(Int, n)), agents
end

"modifies `agents`"
function step!(A, nodes, agents, p, τ, τr, m, q, g)
	infection_count = 0

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
					infection_count += 1
				end
			end

			recovered_neighbors = [b for b in agents if (b.location == a.location) && (b.state == recovered)]
			for b in recovered_neighbors
				if τr > rand()
					b.state = infectious
					b.days_infected = 0
					infection_count += 1
				end
			end
		end
	end

	# check quarntines
	for (i, n) in enumerate(nodes)
		if n.quarantine_state == notquarantined
			node_infectious_agents = count(a -> (a.state == infectious) && (a.location == i), agents)
			node_total_agents = count(a -> a.location == i, agents)

			if (node_infectious_agents / node_total_agents) > g
				n.quarantine_state = quarantined
			end
		else
			if n.days_quarantined > q
				n.quarantine_state = notquarantined
				n.days_quarantined = 0
			end
		end
	end

	# travel
	for a in agents
		if (m > rand()) && (nodes[a.location].quarantine_state == notquarantined)
			# agent travels
			potential_destinations = findall(==(1), A[a.location, :])
			destination_weights = [nodes[i].degree for i in potential_destinations if nodes[i].quarantine_state == notquarantined]
			if !isempty(destination_weights)
				destination = sample(potential_destinations, Weights(destination_weights))
				a.location = destination
			end
		end
	end

	# update days infected
	for a in agents
		if a.state == infectious
			a.days_infected += 1
		end
	end

	# update quarntines
	for n in nodes
		if n.quarantine_state == quarantined
			n.days_quarantined += 1
		end
	end

	return infection_count
end

function infectiouscount(agents)
	count(a -> a.state == infectious, agents)
end

#=
trials = 100
steps = 200
d = 0.2
p = 7
τ = 0.3
τr = τ/10
m = 0.1
q = 28
g = 0.25
Is = []
for i in 1:trials
	i |> println

	A, nodes, agents = init(100, d, 10)

	I = zeros(Int, steps)
	I[1] = 1
	for j in 2:steps
		step!(A, nodes, agents, p, τ, τr, m, q, g)
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
d = 0.2
p = 7
τ = 0.15
τr = τ/10
m = 0.05
q = 20
g = 0.7
Is = []
for i in 1:trials
	i |> println

	A, nodes, agents = init(100, d, 10)
	
	I = zeros(Int, steps)
	I[1] = 1
	for j in 2:steps
		step!(A, nodes, agents, p, τ, τr, m, q, g)
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
=#

ds = range(0.1, 0.3; length=10)  # connectedness of network
ps = range(5, 14; step=1)        # infectious days
τs = range(0.05, 0.3; length=10) # infectious rate
τrs = τs ./ 3                    # reinfectious rate
ms = range(0.05, 0.2; length=10) # transportation rate
qs = ps .* 2                     # quarantine length
gs = range(0.1, 0.5; length=10)  # quarantine threshold

default_d = 0.2
default_p = 10
default_τ = 0.2
default_τr = default_τ / 3
default_m = 0.1
default_q = 20
default_g = 0.25

steps = 100

# (d, p)
I = zeros(Int, 10, 10)
for (i, d) in enumerate(ds)
	for (j, p) in enumerate(ps)
		A, nodes, agents = init(100, d, 10)

		total_infections = 1
		for _ in 2:steps
			total_infections += step!(A, nodes, agents, p, default_τ, default_τr, default_m, default_q, default_g)
		end

		I[i, j] = total_infections
	end
end

p1 = heatmap(ps, ds, I, xlabel="\$ p\$", ylabel="\$ d\$")

# (q, g)
I = zeros(Int, 10, 10)
for (i, q) in enumerate(qs)
	for (j, g) in enumerate(gs)
		A, nodes, agents = init(100, default_d, 10)

		total_infections = 1
		for _ in 2:steps
			total_infections += step!(A, nodes, agents, default_p, default_τ, default_τr, default_m, q, g)
		end

		I[i, j] = total_infections
	end
end

p2 = heatmap(gs, qs, I, xlabel="\$ g\$", ylabel="\$ q\$")

# (τ, τr)
I = zeros(Int, 10, 10)
for (i, τ) in enumerate(τs)
	for (j, τr) in enumerate(τrs)
		A, nodes, agents = init(100, default_d, 10)

		total_infections = 1
		for _ in 2:steps
			total_infections += step!(A, nodes, agents, default_p, τ, τr, default_m, default_q, default_g)
		end

		I[i, j] = total_infections
	end
end

p3 = heatmap(τrs, τs, I, xlabel="\$ τ_r\$", ylabel="\$ τ\$")

# (p, q)
I = zeros(Int, 10, 10)
for (i, p) in enumerate(ps)
	for (j, q) in enumerate(qs)
		A, nodes, agents = init(100, default_d, 10)

		total_infections = 1
		for _ in 2:steps
			total_infections += step!(A, nodes, agents, p, default_τ, default_τr, default_m, q, default_g)
		end

		I[i, j] = total_infections
	end
end

p4 = heatmap(qs, ps, I, xlabel="\$ q\$", ylabel="\$ p\$")

plot(p1, p2, p3, p4)
savefig("plots/10.14.png")
