using ITensors
using ITensorMPS

#include("cliffordaug.jl")
#using .cliffordaug
begin
    


N = 5                  # Number of spins
J = 1.0                 # Heisenberg coupling
dt = 0.1               # Time step
tmax = 2              # Total evolution time

maxdim = 100            # Maximum MPS bond dimension
cutoff = 1E-10          # SVD truncation cutoff

sites=siteinds("S=1/2",5)
psi=random_mps(sites)
entanglement(psi,2)
    
h=OpSum()

for j in 1:N-1
    global h += J, "Sz", j, "Sz", j+1
    global h += J/2, "S+", j, "S-", j+1
    global h += J/2, "S-", j, "S+", j+1
end

H = MPO(h, sites)

times = Float64[]
entropies = Float64[]

middle = div(N, 2)
nsteps=floor(tmax/dt)  
for step in 1:nsteps

    # Real-time TDVP evolution
    global psi = tdvp(
        H,
        -1im * dt,
        psi;
        nsteps = 1,
        maxdim = maxdim,
        cutoff = cutoff,
        nsite = 2,
        normalize=true
    )

    # Measure entanglement
    entropy = entanglement(psi, middle)

    push!(times, step*dt)
    push!(entropies, entropy)

    println("t = $(round(step*dt,digits=3)), S = $entropy")
end
methods(cliffordaug.optimization)
methods(cliffordaug.cliffordaug_sites)
beta=optimization(psi)
end