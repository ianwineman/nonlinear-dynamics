using FFTW, Statistics

function logistic_map(x; r=3.25)#r=3.58)
	return r*x*(1-x)
end

H(probs) = -sum([p*log(p) for p in probs])

# Periodic timeseries
r = 3.25
u0 = 0.4
traj = [u0]

for _ in 1:1000-1
	push!(traj, logistic_map(traj[end]; r=r))
end

probs = abs2.(fft(traj .- mean(traj))) 
probs ./= sum(probs)

printstyled("Periodic timeseries spectral entropy: "; color=:light_magenta)
printstyled("$(round(H(probs); digits=4))\n"; color=:light_cyan)

# Chaotic timeseries
r = 3.58
u0 = 0.4
traj = [u0]

for _ in 1:1000-1
	push!(traj, logistic_map(traj[end]; r=r))
end

probs = abs2.(fft(traj .- mean(traj))) 
probs ./= sum(probs)

printstyled("Chaotic  timeseries spectral entropy: "; color=:light_magenta)
printstyled("$(round(H(probs); digits=4))\n"; color=:light_cyan)

# Noise timeseries
traj = rand(1000)

probs = abs2.(fft(traj .- mean(traj))) 
probs ./= sum(probs)

printstyled("Pure noise          spectral entropy: "; color=:light_magenta)
printstyled("$(round(H(probs); digits=4))\n"; color=:light_cyan)

# Theoretical maximum entropy
printstyled("Theoretical maximum spectral entropy: "; color=:light_magenta)
printstyled("$(round(log(1000); digits=4))\n"; color=:light_cyan)