using Plots

function logistic_map(x)
	return 4*x*(1-x)
end

function tent_map(x)
	if x ≤ 0.5
		return 2*x
	else
		return 2 - 2*x
	end
end

function transformed_logistic_map(x)
	return (2/π)*asin(sqrt(logistic_map(sin((π*x)/2)^2)))
end

function trajectory(f, x0, N)
	t = zeros(typeof(x0), N)
	t[1] = x0
	for n in 2:N
		t[n] = f(t[n-1])
	end
	return t
end

x0 = π/4
N = 50

t = trajectory(tent_map, x0, N)
plot(1:N, t, 
	xlabel="\$n\$", ylabel="\$x_n\$",
	xlim=[0,N], ylim=[0,1],
	label="Tent map: \$μ=2\$",
	lc=:blue
)

t = trajectory(transformed_logistic_map, x0, N)
plot!(1:N, t, 
	xlabel="\$n\$", ylabel="\$x_n\$",
	xlim=[0,N], ylim=[0,1],
	label="Transformed logistic map: \$r=4\$",
	ls=:dash, lc=:orange
)
savefig("plots/exercise5_4.png")
