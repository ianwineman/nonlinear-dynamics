using Plots, FFTW, Statistics, DifferentialEquations, DSP

function brusselator!(du, u, p, t)
   a, b = p
   x, y = u
   du[1] = 1 - (b+1)*x + a*(x^2)*y
   du[2] = b*x - a*(x^2)*y
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

prob = ODEProblem(brusselator!, [0.5909924611496271,4.843385159662832], (0.0,13.0), [0.3, 1.5])
sol = solve(prob, dtmax=0.01)

u = sol[1, :] .- mean(sol[1, :])
v = hilbert(u)
z = u .+ v*1im
z = z ./ (abs.(z))
ψ = atan.(imag.(z) ./ real.(z)) .+ shift.(real.(z), imag.(z))

p1 = plot(sol, idxs=(1,2), xlabel="\$ u\$", ylabel="\$ w\$", label=false, aspect_ratio=1, title="Brusselator \$ (u(t),w(t))\$")
p2 = plot(sol[1, :], xlabel="\$ t\$", ylabel="\$ u\$", label=false, title="\$ u(t)\$")
p3 = plot(z, xlabel="\$ Re(z/|z|)\$", ylabel="\$ Im(z/|z|)\$", label="\$ z(t)=u(t) + iH[u(t)]\$", aspect_ratio=1, title="Analytic signal \$ z(t)\$")
p4 = plot(ψ, xlabel="\$ t\$", ylabel="\$ ψ\$", label=false, title="Phase angle \$ ψ(t)\$")
plot(p1,p2,p3,p4,layout=(2,2),size=(800,600), linewidth=2, linecolor=[Colors.JULIA_LOGO_COLORS.blue Colors.JULIA_LOGO_COLORS.red Colors.JULIA_LOGO_COLORS.green Colors.JULIA_LOGO_COLORS.purple])
savefig("plots/exercise2_9.png")
