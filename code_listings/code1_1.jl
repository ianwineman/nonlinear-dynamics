# Code 1.1 (Datseris & Parlitz)
using DynamicalSystems

function lorenz_rule(u, p, t)
    σ, ρ, β = p 
    x, y, z = u 
    dx = σ*(y-x)
    dy = x*(ρ-z)-y
    dz = x*y - β*z
    return SVector(dx,dy,dz) # construct a statically-sized vector SVector
end

p  = [10.0, 28.0, 8/3]
uθ = [ 0.0, 10.0, 0.0]
lorenz = ContinuousDynamicalSystem(lorenz_rule, uθ, p)

T = 100.0 # total time
Δt = 0.01 # sampling time
A = trajectory(lorenz, T; Δt)