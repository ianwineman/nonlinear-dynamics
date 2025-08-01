using Plots, DifferentialEquations

function popdy!(du, u, p, t)
   r1, w1, K1, r2, w2, K2 = p
   n1, n2 = u
   du[1] = r1*(1-n1)*n1 + w1*n1*n2*K2
   du[2] = r2*(1-n2)*n2 + w2*n1*n2*K1
end

function equilibria(r1, w1, K1, r2, w2, K2)
	n1_star = ((r1*r2) + (r2*w1*K2))/((r1*r2) - (w1*w2*K1*K2))
	n2_star = ((r1*r2) + (r1*w2*K1))/((r1*r2) - (w1*w2*K1*K2))

	println((n1_star, n2_star))
	return (n1_star, n2_star)
end

prob = ODEProblem(popdy!, [0.000000001,0.0], (0.0,50.0), [1.5,-0.1,10.0,1.2,-0.1,10.0])
sol = solve(prob, dtmax=0.01)

equilibria([1.5,-0.1,10.0,1.2,-0.1,10.0]...)
plot(sol, ylim=[0,1])
