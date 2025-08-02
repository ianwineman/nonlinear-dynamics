using Plots, FFTW, Statistics, DifferentialEquations, DSP, LinearAlgebra

function lorenz63!(du, u, p, t)
	σ, ρ, β = p
    x, y, z = u

    du[1] = σ * (y - x)
    du[2] = -x * z + ρ * x - y
    du[3] = x * y - β * z
end

function shift(x,y)
    if x>=0 && y>=0
        return 0
    elseif (x<0 && y>=0) || (x<0 && y<0)
        return π
    else
        return 2π
    end
end

function phase(sol) 
    period = 3.46
    t = sol.t
    ϕ = ((t ./ period) .* 2π) .% 2π
end

function angle(sol) 
    x, y = sol[1, :], sol[2, :]
    θ = atan.(y ./ x) .+ shift.(x, y)
end

function arclength(sol)
    steps = [(i==length(sol.u)) ? abs(norm(sol.u[1] .- sol.u[i])) : abs(norm(sol.u[i+1] .- sol.u[i]))  for (i,_) in enumerate(sol.u)]
    arclengths = [sum(steps[1:i]) for i=1:length(sol.u)]
    normalized_arclengths = arclengths ./ arclengths[end]
end

function hilberphaseangle(sol)
    u = sol[1, :] .- mean(sol[1, :])
    v = hilbert(u)
    z = u .+ v*1im
    z = z ./ (abs.(z))
    ψ = atan.(imag.(z) ./ real.(z)) .+ shift.(real.(z), imag.(z))
end

prob = ODEProblem(lorenz63!, [-9.956186631140996,-27.19316994087894,104.00503095666596], (0.0,3.46), [10.0, 160.0, 8/3])
sol = solve(prob, dtmax=0.01)

p1 = plot(sol, idxs=(1,2), xlabel="\$ x\$", ylabel="\$ y\$", label=false, title="Lorenz-63")
p2 = plot(phase(sol), xlabel="\$ t\$", ylabel="\$ ϕ\$", label=false, title="Phase", xlim=[0,345], ytick=([0.0,π,2π], ["0", "π", "2π"]))
p3 = plot(angle(sol), xlabel="\$ t\$", ylabel="\$ θ\$", label=false, title="Angle about \$ (0,0)\$")
p4 = plot(arclength(sol), xlabel="\$ t\$", ylabel="\$ l\$", label=false, title="Arclength")
p5 = plot(hilberphaseangle(sol), xlabel="\$ t\$", ylabel="\$ ψ\$", label=false, title="Phase angle for \$ x(t)\$")
plot(p1,p2,p3,p4,p5, size=(800,600), linewidth=2, linecolor=[:lightskyblue Colors.JULIA_LOGO_COLORS.blue Colors.JULIA_LOGO_COLORS.red Colors.JULIA_LOGO_COLORS.green Colors.JULIA_LOGO_COLORS.purple], layout=@layout [a; grid(2,2)])
savefig("plots/exercise2_10.png")
