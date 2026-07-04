using DifferentialEquations
using Statistics
using LinearAlgebra
using Plots

function sprott!(du, u, p, t)
   x, y, z = u
   du[1] = y + 2*x*y + x*z
   du[2] = 1 - 2*x^2 + y*z
   du[3] = x - x^2 - y^2
end

"u = solv.u; dim(u) = 3 only"
function density(u, ϵ)
   r  = u .|> norm |> maximum
   bc = r / ϵ |> round |> Int # box count from origin in each direction

   xrange = range(-r, r; length=bc*2)
   yrange = range(-r, r; length=bc*2)
   zrange = range(-r, r; length=bc*2)

   boxes = [[x, y, z] for x in xrange for y in yrange for z in zrange]

   density = zeros(boxes |> length)
   for v in u
      box_dis = [norm(b - v) for b in boxes]
      i = argmin(box_dis)
      density[i] += 1
   end
   density ./= length(u)
   return boxes, density
end

prob = ODEProblem(sprott!, [1.0, 0.0, 0.0], (0.0, 500.0), [])
solv = solve(prob; dtmax=0.1)

b, d = density(solv.u, 0.1)

w1 = 0 # space average
w2 = (getindex.(solv.u, 2) + getindex.(solv.u, 3)) |> mean # time average
for (box, box_dens) = zip(b, d)
   y, z = box[2:3]
   global w1 += ((y + z) * box_dens)
end

w1 = round(w1; digits=5)
w2 = round(w2; digits=5)

p1 = plot(
   solv, 
   idxs=(1, 2, 3),
   label="\$ (y+z) = $w1\$ (space average)\n\$ (y+z) = $w2\$ (time average)",
   legendtitle="Average volume dissipation rate",
   legendposition=:outerbottom,
   title="Sprott \$ u_0=(1, 0, 0)\$",
   #lα=0.25,
   axis=false,
   grid=false
)

prob = ODEProblem(sprott!, [2.0, 0.0, 0.0], (0.0, 500.0), [])
solv = solve(prob; dtmax=0.1)

b, d = density(solv.u, 0.1)

w1 = 0 # space average
w2 = (getindex.(solv.u, 2) + getindex.(solv.u, 3)) |> mean # time average
for (box, box_dens) = zip(b, d)
   y, z = box[2:3]
   global w1 += ((y + z) * box_dens)
end

w1 = round(w1; digits=5)
w2 = round(w2; digits=5)

p2 = plot(
   solv, 
   idxs=(1, 2, 3),
   label="\$ (y+z) = $w1\$ (space average)\n\$ (y+z) = $w2\$ (time average)",
   legendtitle="Average volume dissipation rate",
   legendposition=:outerbottom,
   title="Sprott \$ u_0=(2, 0, 0)\$",
   #lα=0.25,
   lc=:purple,
   axis=false,
   grid=false
)

plot(p1, p2, size=(800, 400))
savefig("plots/8.16.png")
