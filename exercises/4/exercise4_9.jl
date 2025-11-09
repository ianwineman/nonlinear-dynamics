using Plots, LinearAlgebra

function h(z)
	return z^3 - 1
end

h_J_inv = z -> 1/(3*z^2)

function damped_newtons_method(f, X0, J_inv; atol=0.001)
	Xj = X0
	while norm(f(Xj)) > atol
		δ = 1.0
		f_Xj = f(Xj)

		while norm(f(Xj .- δ .* J_inv(Xj)*f_Xj)) >= norm(f_Xj)
			δ /= 2
			if δ < atol
				break
			end
		end

		next_Xj = Xj .- δ .* J_inv(Xj) * f_Xj
		Xj = next_Xj
	end
	return Xj
end

fp1 = 1 + 0im
fp2 = -1/2 + √3/2im
fp3 = -1/2 - √3/2im
x0s = -1.5:0.01:1.5
y0s = -1.5:0.01:1.5

fp1_basin = []
fp2_basin = []
fp3_basin = []

for x in x0s
	for y in y0s
		z = complex(x,y)
		dnm = damped_newtons_method(h, z, h_J_inv)
		if abs(fp1 - dnm) < 0.001
			push!(fp1_basin, (x,y))
		elseif abs(fp2 - dnm) < 0.001
			push!(fp2_basin, (x,y))
		elseif abs(fp3 - dnm) < 0.001
			push!(fp3_basin, (x,y))
		end
	end
end

scatter(first.(fp1_basin), last.(fp1_basin), msw=0, ms=1, mc=:red, label=false, shape=:rect)
scatter!(first.(fp2_basin), last.(fp2_basin), msw=0, ms=1, mc=:blue, label=false, shape=:rect)
scatter!(first.(fp3_basin), last.(fp3_basin), msw=0, ms=1, mc=:orange, label=false, shape=:rect)

plot!(
	title="Newton fractal: \$h(z) = z^3 -1\$", 
	grid=false, 
	xlim=[-1.5,1.5],
	ylim=[-1.5,1.5],
	xtick=[-1.5,0.0,1.5],
	ytick=[-1.5,0.0,1.5],
)
savefig("plots/exercise4_9.png")
