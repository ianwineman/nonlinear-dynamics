using Plots

function model_f(x; p=[0.0, 1.0])
	ϵ, ρ = p
	ϵ + ρ*x - x^3
end

function model_J(x; p=[0.0, 1.0])
	ϵ, ρ = p
	ρ - 3*x^2
end

function model_J_inv(x; p=[0.0, 1.0])
	ϵ, ρ = p
	1/(ρ - 3*x^2)
end

function damped_newtons_method(f, x0, J_inv; p=[], atol=0.001)
	xj = x0
	while abs(f(xj; p=p)) > atol
		δ = 1.0
		f_xj = f(xj; p=p)

		while abs(f(xj - δ*J_inv(xj; p=p)*f_xj; p=p)) >= abs(f_xj)
			δ /= 2
			if δ < atol
				break
			end
		end

		next_xj = xj - δ*J_inv(xj; p=p)*f_xj
		xj = next_xj
	end
	return xj
end

function continuation(f, zm1, zm0, N, J_inv; p=[])
	# z = [p,x]
	zm, zm_min1 = zm1, zm0
	zms = [zm, zm_min1]
	for m in 1:N
		ẑ = zm .+ (zm .- zm_min1)                                                          # predictor (secant)
		zm_pls1 = [first(ẑ), damped_newtons_method(f, last(ẑ), J_inv; p=[p[1], first(ẑ)])] # corrector (Newton's method)
		push!(zms, zm_pls1)
		zm_min1 = zm
		zm = zm_pls1
	end
	return zms
end

let
	anim = @animate for ϵ in vcat([0.0 for _ in 1:9], collect(0.0:0.01:0.5), collect(0.49:-0.01:0.0), [0.0 for _ in 1:8], collect(0.0:-0.01:-0.5))
		println("ϵ=$ϵ")
		stable_fixed_points = []
		unstable_fixed_points = []

		# continuation 1
		c1_z0 = [5.00,damped_newtons_method(model_f, 5.0, model_J_inv; p=[ϵ, 5.00])]
		c1_z1 = [4.99,damped_newtons_method(model_f, 5.0, model_J_inv; p=[ϵ, 4.99])]
		c1 = continuation(model_f, c1_z1, c1_z0, 1000, model_J_inv; p=[ϵ,c1_z0[1]])

		# continuation 2
		c2_z0 = [5.00,damped_newtons_method(model_f, -5.0, model_J_inv; p=[ϵ, 5.00])]
		c2_z1 = [4.99,damped_newtons_method(model_f, -5.0, model_J_inv; p=[ϵ, 4.99])]
		c2 = continuation(model_f, c2_z1, c2_z0, 1000, model_J_inv; p=[ϵ,c2_z0[1]])

		# continuation 3
		c3_z0 = [5.00,damped_newtons_method(model_f, 0.0, model_J_inv; p=[ϵ, 5.00])]
		c3_z1 = [4.99,damped_newtons_method(model_f, 0.0, model_J_inv; p=[ϵ, 4.99])]
		c3 = continuation(model_f, c3_z1, c3_z0, 1000, model_J_inv; p=[ϵ,c3_z0[1]])

		for fixed_point in vcat(c1,c2,c3)
			ρ, x = fixed_point
			if model_J(x; p=[ϵ,ρ]) < 0.0
				push!(stable_fixed_points, (ρ, x))
			else
				push!(unstable_fixed_points, (ρ, x))
			end
		end
		
		scatter(
			first.(stable_fixed_points), last.(stable_fixed_points), 
			xlabel="\$ρ\$", ylabel="\$x^*\$",
			legendtitle="ϵ=$ϵ",
			label="Stable",
			msw=0, ms=3, mc=:blue,
			title="Bifurcation Diagram: \$ẋ = ϵ + ρx - x^3\$"
		)
		scatter!(
			first.(unstable_fixed_points), last.(unstable_fixed_points), 
			label="Unstable",
			msw=0, ms=3, mc=:red
		)
	end
	gif(anim, "plots/exercise4_3.gif", fps=10)
end
