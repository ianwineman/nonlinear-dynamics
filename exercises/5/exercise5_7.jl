using LinearAlgebra, DifferentialEquations, Plots

function correlation_sum(X, ϵ, w)
	N = length(X)
	C = 0

	for i in 1:(N - w - 1)
		for j in (1 + w + i):N
			C += ifelse(norm(X[i] .- X[j]) < ϵ, 1, 0)
		end
	end

	return (2*C)/((N - w)*(N - w - 1))
end

function standard_map(θ, v; k=1.0)
	return mod.([θ + v + k*sin(θ), v + k*sin(θ)], 2π)
end

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function henonheiles!(du, u, p, t)
   x, y, vx, vy = u
   du[1] = vx
   du[2] = vy
   du[3] = -x - 2*x*y
   du[4] = -y - x^2 + y^2
end


traj = [[0.1,0.1]]
for _ = 1:1_000
	push!(traj, standard_map(traj[end]...))
end

ϵs = [exp(i) for i in -10.0:0.1:3.0]
Cs = [correlation_sum(traj, ϵ, 3) for ϵ in ϵs]

n = findfirst(x->!isinf(x) && !isnan(x), log.(Cs))
plot(
	log.(ϵs)[n:end], 
	log.(Cs)[n:end],
	xlabel="\$\\log{(ϵ)}\$",
	xtick=-6:1:3,
	ylabel="\$\\log{(C)}\$",
	title="Standard Map",
	label=false,
	lw=2,
	lc=:blue
)

xstart, xend = -4.0, -1.0
dim = (log.(Cs)[findfirst(x->x==xend, -10.0:0.1:3.0)] - log.(Cs)[findfirst(x->x==xstart, -10.0:0.1:3.0)])/(xend - xstart)
plot!(
	range(xstart,xend; length=100),
	range(xstart,xend; length=100) .* dim .+ log.(Cs)[findfirst(x->x==0.0, -10.0:0.1:3.0)],
	label="\$Δ^{(C)} ≈ $(round(dim; digits=2))\$",
	lc=:purple,
	lw=2
)
savefig("plots/exercise5_7_sm.png")


p = [10.0, 28.0, 8/3]
u0 = [20.0, 20.0, 20.0]
prob = ODEProblem(lorenz63!, u0, (0.0,200.0), p)
solv = solve(prob)

ϵs = [exp(i) for i in -20.0:0.1:5.0]
Cs = [correlation_sum(solv.u, ϵ, 3) for ϵ in ϵs]

n = findfirst(x->!isinf(x) && !isnan(x), log.(Cs))
plot(
	log.(ϵs)[n:end], 
	log.(Cs)[n:end],
	xlabel="\$\\log{(ϵ)}\$",
	xtick=-4:1:5,
	ylabel="\$\\log{(C)}\$",
	title="Lorenz-63",
	label=false,
	lw=2,
	lc=:blue
)

xstart, xend = -3.5, 2.0
dim = (log.(Cs)[findfirst(x->x==xend, -20.0:0.1:5.0)] - log.(Cs)[findfirst(x->x==xstart, -20.0:0.1:5.0)])/(xend - xstart)
plot!(
	range(xstart,xend; length=100),
	range(xstart,xend; length=100) .* dim .+ log.(Cs)[findfirst(x->x==0.0, -20.0:0.1:5.0)] .- 0.5,
	label="\$Δ^{(C)} ≈ $(round(dim; digits=2))\$",
	lc=:purple,
	lw=2
)
savefig("plots/exercise5_7_lo.png")


u0 = [0.0, -0.31, 0.354, 0.059]
prob = ODEProblem(henonheiles!, u0, (0.0,300.0), [])
solv = solve(prob)

ϵs = [exp(i) for i in -20.0:0.01:1.0]
Cs = [correlation_sum(solv.u, ϵ, 3) for ϵ in ϵs]

n = findfirst(x->!isinf(x) && !isnan(x), log.(Cs))
dim = correlation_dimension(solv.u, exp(-1), 3)
plot(
	log.(ϵs)[n:end], 
	log.(Cs)[n:end],
	xlabel="\$\\log{(ϵ)}\$",
	ylabel="\$\\log{(C)}\$",
	title="Hénon-Heiles",
	label=false,
	lw=2,
	lc=:blue
)

xstart, xend = -3.0, -0.5
dim = (log.(Cs)[findfirst(x->x==xend, -20.0:0.01:1.0)] - log.(Cs)[findfirst(x->x==xstart, -20.0:0.01:1.0)])/(xend - xstart)
plot!(
	range(xstart,xend; length=100),
	range(xstart,xend; length=100) .* dim .+ log.(Cs)[findfirst(x->x==0.0, -20.0:0.01:1.0)] .- 0.25,
	label="\$Δ^{(C)} ≈ $(round(dim; digits=2))\$",
	lc=:purple,
	lw=2
)
savefig("plots/exercise5_7_hh.png")
