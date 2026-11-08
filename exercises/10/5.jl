using DifferentialEquations, StaticArrays

### from book ###
function Roessler_ODE(u,p,t)
    a, b, c = p
    x, y, z = u
    dx =  - y - z
    dy =  x + a*y
    dz = b + z*(x-c)
    return SVector(dx, dy, dz)
end

function Roessler_MSF_ODE(u,p,t)
    a, b, c, alpha, Jh = p
    x1, x2, x3, y1, y2, y3  = u
    dx1 =  - x2 - x3
    dx2 =  x1 + a*x2
    dx3 = b + x3*(x1-c)

   # linearized equation
    Jg = [0 -1 -1 ; 1 a 0 ; x3 0 x1-c]
    dy = (Jg .- alpha.*Jh) * u[4:6]

    return SVector(dx1, dx2, dx3, dy[1], dy[2], dy[3])
end

function MSF(uinit,a,b,c,Jh)
  
    nle = 20000 # number of renormalization steps
    tle = 6.  # period of time between renormalizations

    alpha_min = -2.
    alpha_max = 8.
    nalpha = 200
    alpha_vec = LinRange(alpha_min, alpha_max, nalpha)
  #  avec = zeros(nalpha)   # alpha values
    LE_vec = zeros(nalpha)  # corresponding Lyapunov exponent
    
   # uinit = u0[1:3]
    
    for ialpha = 1:nalpha
    
        alpha = alpha_vec[ialpha] 
        cum = 0.

        p = a, b, c, alpha, Jh     
        u0 = [uinit ; SVector(1., 0., 0.) ]           

        for ile = 1:nle   # loop over renormalization steps

            tspan = (0., tle)
            prob = ODEProblem(Roessler_MSF_ODE, u0, tspan, p)
            sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8,tsamp = tle)
    
            yvec = [u[4:6] for u in sol.u]  # perturbation vector
            yend = yvec[end]
            ynorm = sqrt( sum(yend .* yend)  )
            yend = yend / ynorm      # renormalization
            cum = cum + log(ynorm)  
    
            uvec = [u[1:3] for u in sol.u]
            uend = uvec[end]
    
            u0 = [uend ; yend]   # initial condition for next integration step
    
        end
    
        lambda_max = cum / (nle*tle) # Lyapunov exponent
       # println(alpha, "  ", lambda_max)
    
       # avec[ialpha] = alpha
        LE_vec[ialpha] = lambda_max
    
    end # of alpha loop
    return alpha_vec, LE_vec
end # of function

function MSF(uinit,a,b,c,Jh)
  
    nle = 20000 # number of renormalization steps
    tle = 6.  # period of time between renormalizations

    alpha_min = -2.
    alpha_max = 8.
    nalpha = 200
    alpha_vec = LinRange(alpha_min, alpha_max, nalpha)
  #  avec = zeros(nalpha)   # alpha values
    LE_vec = zeros(nalpha)  # corresponding Lyapunov exponent
    
   # uinit = u0[1:3]
    
    for ialpha = 1:nalpha
    
        alpha = alpha_vec[ialpha] 
        cum = 0.

        p = a, b, c, alpha, Jh     
        #u0 = [uinit ; SVector(1., 0., 0.) ]
        u0 = SVector{6}([uinit ; SVector(1., 0., 0.)])           

        for ile = 1:nle   # loop over renormalization steps

            tspan = (0., tle)
            prob = ODEProblem(Roessler_MSF_ODE, u0, tspan, p)
            sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
            yvec = [u[4:6] for u in sol.u]  # perturbation vector
            yend = yvec[end]
            ynorm = sqrt( sum(yend .* yend)  )
            yend = yend / ynorm      # renormalization
            cum = cum + log(ynorm)  
    
            uvec = [u[1:3] for u in sol.u]
            uend = uvec[end]
    
            #u0 = [uend ; yend]   # initial condition for next integration step
    		u0 = SVector{6}([uend ; yend])
        end
    
        lambda_max = cum / (nle*tle) # Lyapunov exponent
       # println(alpha, "  ", lambda_max)
    
       # avec[ialpha] = alpha
        LE_vec[ialpha] = lambda_max
    
    end # of alpha loop
    return alpha_vec, LE_vec
end # of function

a = 0.2
b = 0.2 
c = 5.7 

p = a, b, c
u0 = SVector(0.1, 0.1, 0.0)

tspan = (0.0, 1000)
prob = ODEProblem(Roessler_ODE, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

u = sol.u
uinit = u[end]

Jh = [1 0 0; 0 0 0; 0 0 0] # x drives x-ODE
alpha_vec, LE_vec_xx = MSF(uinit, a, b, c, Jh)
### end from book ###

using LinearAlgebra
using Plots

# ring
W = [0 1 0 0 0 0 0 0 0 1;
     1 0 1 0 0 0 0 0 0 0;
     0 1 0 1 0 0 0 0 0 0;
     0 0 1 0 1 0 0 0 0 0;
     0 0 0 1 0 1 0 0 0 0;
     0 0 0 0 1 0 1 0 0 0;
     0 0 0 0 0 1 0 1 0 0;
     0 0 0 0 0 0 1 0 1 0;
     0 0 0 0 0 0 0 1 0 1;
     1 0 0 0 0 0 0 0 1 0;
]

G = zeros(10, 10)

for i in 1:10
	for j in 1:10
		if i == j
			G[i, j] = sum(W[i, :])
		else
			G[i, j] = - W[i, j]
		end
	end
end

λ_ring = eigen(G).values

# alpha_vec[127] = 4.331658291457286
σmax_ring = alpha_vec[127] / λ_ring[end]

plot(
	alpha_vec, 
	LE_vec_xx, 
	markershape=:circ, 
	markersize=1, 
	label=false, 
	xlabel="\$ α\$", 
	ylabel="\$ λ_1\$"
)
vline!(
	[alpha_vec[44]], # α = 0.161
	label="\$α = $(round(alpha_vec[44]; digits=3))\$"
)
vline!(
	[alpha_vec[127]], # α = 4.332
	label="\$ α = $(round(alpha_vec[127]; digits=3))\$"
)
scatter!(λ_ring, zeros(length(λ_ring)), label="\$ σ=1\$ (ring), sync when \$ σ \\ ∈ \\ (0, $(round(σmax_ring; digits=3)))\$")

# starred ring
W = [0 1 0 0 0 0 0 0 0 1 1;
     1 0 1 0 0 0 0 0 0 0 1;
     0 1 0 1 0 0 0 0 0 0 1;
     0 0 1 0 1 0 0 0 0 0 1;
     0 0 0 1 0 1 0 0 0 0 1;
     0 0 0 0 1 0 1 0 0 0 1;
     0 0 0 0 0 1 0 1 0 0 1;
     0 0 0 0 0 0 1 0 1 0 1;
     0 0 0 0 0 0 0 1 0 1 1;
     1 0 0 0 0 0 0 0 1 0 1;
     1 1 1 1 1 1 1 1 1 1 0;
]

G = zeros(11, 11)

for i in 1:11
	for j in 1:11
		if i == j
			G[i, j] = sum(W[i, :])
		else
			G[i, j] = - W[i, j]
		end
	end
end

λ_starredring = eigen(G).values

# alpha_vec[127] = 4.331658291457286
σmax_starredring = alpha_vec[127] / λ_starredring[end]

scatter!(0.25 .* λ_starredring, zeros(length(λ_starredring)), label="\$ σ=0.25\$ (starred ring), sync when \$ σ \\ ∈ \\ (0, $(round(σmax_starredring; digits=3)))\$")
savefig("plots/10.5.png")
