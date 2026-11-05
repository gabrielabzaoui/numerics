using ITensors
using ITensorMPS

begin
include("cliffordaug.jl")
using .cliffordaug
    


N = 16                 # Number of spins
J = 1.0                 # Heisenberg coupling
dt = 0.1               # Time step
tmax = 4              # Total evolution time

maxdim = 100            # Maximum MPS bond dimension
cutoff = 1E-10          # SVD truncation cutoff

sites=siteinds("S=1/2",N)
psi=random_mps(sites)
cliffordaug.entanglement(psi,2)
    
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
    #psi_opt = cliffordaug.optimization(psi,4)

    #println(
    #"S diff = $(cliffordaug.entanglement(psi, middle)-cliffordaug.entanglement(psi_opt, middle))"
    #)
   println(step)
end
S_before = cliffordaug.entanglement(psi,middle)
    psi_opt = cliffordaug.optimization(deepcopy(psi), 2)
    S_after = cliffordaug.entanglement(psi_opt,middle)

    println("Before: $S_before")
    println("After:  $S_after")
    println("Diff:   $(S_before - S_after)")
    

end