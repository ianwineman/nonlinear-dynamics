using LinearAlgebra
using Plots
using DifferentialEquations
using ProgressMeter

# x → dx/dt
function coupledrossler!(du, u, p, t)
	a, b, c, W, σ = p

	x = u[:, 1]
	y = u[:, 2]
	z = u[:, 3]

	N = first(size(W))
    for i in 1:N
    	du[i, 1] = - u[i, 2] - u[i, 3] + σ * sum([W[i, j] * (u[j, 1] - u[i, 1]) for j in 1:N])
    	du[i, 2] = u[i, 1] + (a * u[i, 2])
    	du[i, 3] = b + u[i, 3] * (u[i, 1] - c)
    end
end

W = [0 7 0 5 0;
     7 0 3 0 4;
     0 3 0 0 6;
     5 0 0 0 0;
     0 4 6 0 0]

p = (0.2, 0.2, 5.7, W, 0.15)
u0 = rand(first(size(W)), 3) .* 10.0

prob = ODEProblem(coupledrossler!, u0, (0.0, 50.0), p)
solv = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-10)

x1 = solv[1, 1, :]
y1 = solv[1, 2, :]
z1 = solv[1, 3, :]

p1 = plot(
	x1, 
	y1, 
	z1, 
	label=false,
	xlabel="\$ x\$",
	ylabel="\$ y\$",
	zlabel="\$ z\$"
)
for i in 2:first(size(W))
	plot!(solv[i, 1, :], solv[i, 2, :], solv[i, 3, :], label=false)
end

p2 = plot(
	solv.t, 
	x1,
	xlabel="\$ t\$",
	ylabel="\$ x\$",
	label=false
)
for i in 2:first(size(W))
	plot!(solv.t, solv[i, 1, :], label=false)
end


# σ sweep
σs = range(0.0, 0.3; length=100)
xerrs = zeros(length(σs))
trials = 100

@showprogress Threads.@threads for (i, σ) in collect(enumerate(σs))
	trial_errs = zeros(trials)

	for trial in 1:trials
		u0 = rand(first(size(W)), 3) .* 10.0
		p = (0.2, 0.2, 5.7, W, σ)
		prob = ODEProblem(coupledrossler!, u0, (0.0, 50.0), p)
		solv = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-10)
		xend = solv[:, 1, :][:, end]
		mn, mx = extrema(xend)
		trial_err = abs(mx - mn)
		trial_errs[trial] = trial_err
	end
	err = sum(trial_errs) / trials
	xerrs[i] = err
end

p3 = plot(
	σs, 
	xerrs,
	xlabel="\$ σ\$",
	ylabel="error",
	ylim=[0, 10],
	ytick=vcat(0:2:10, 1),
	label="\$ x \\ → \\ dx/dt\$",
	lw=2
)


# y → dy/dt
function coupledrossler!(du, u, p, t)
	a, b, c, W, σ = p

	x = u[:, 1]
	y = u[:, 2]
	z = u[:, 3]

	N = first(size(W))
    for i in 1:N
    	du[i, 1] = - u[i, 2] - u[i, 3]
    	du[i, 2] = u[i, 1] + (a * u[i, 2]) + σ * sum([W[i, j] * (u[j, 2] - u[i, 2]) for j in 1:N])
    	du[i, 3] = b + u[i, 3] * (u[i, 1] - c)
    end
end

# σ sweep
σs = range(0.0, 0.3; length=100)
yerrs = zeros(length(σs))
trials = 100

@showprogress Threads.@threads for (i, σ) in collect(enumerate(σs))
	trial_errs = zeros(trials)

	for trial in 1:trials
		u0 = rand(first(size(W)), 3) .* 10.0
		p = (0.2, 0.2, 5.7, W, σ)
		prob = ODEProblem(coupledrossler!, u0, (0.0, 50.0), p)
		solv = solve(prob, Tsit5(); reltol=1e-8, abstol=1e-10)
		yend = solv[:, 2, :][:, end]
		mn, mx = extrema(yend)
		trial_err = abs(mx - mn)
		trial_errs[trial] = trial_err
	end
	err = sum(trial_errs) / trials
	yerrs[i] = err
end

plot!(
	σs, 
	errs,
	xlabel="\$ σ\$",
	ylabel="error",
	label="\$ y \\ → \\ dy/dt\$",
	lw=2
)
hspan!([0.0, 1.0], alpha=0.2, color=:blue, label="sync")

layout = @layout [a b; c]
plot(p1, p2, p3, layout=layout)
savefig("plots/10.4.png")
