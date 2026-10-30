using ITensors
using ITensorMPS
using Revise

let 
    N = 16
    a=2
    J=5
    h=-1
    sites=siteinds("S=1/2",N)
    #site=siteind("S=1/2",N)
    hamiltonian=OpSum()
    ITensors.op(::OpName"Sz",::SiteType"S=1/2") =
    [1 0
     0 -1]
    ITensors.op(::OpName"Sx",::SiteType"S=1/2") =
    [0 1
     1 0]
    
    for i=1:N
        hamiltonian+= -h,"Sx",i
        if i<N
            #Sz=a*op("Sz",sites[i])
            #Sz_prime=a*op("Sz",sites[i+1])
            hamiltonian+= -J,"Sz",i,"Sz",i+1
        
        end   #Sx=a*op("Sx",sites[i])
        
        


        
    end


    nsweeps = 5
    maxdim = [10, 20, 50, 50, 100]
    cutoff = 1e-10
    @time H=MPO(hamiltonian,sites)
    psi0= random_mps(sites;linkdims=10)
    @time energy=dmrg(H,psi0;nsweeps=2,
    maxdim=20,
    cutoff=1e-10)
    @time psi2 = tdvp(H,-1.0im,psi0;nsteps=10)
    println(energy)
    println(inner(psi2',H,psi2))
    #Sz1=op("Sz",sites[1])
    #I2=op("Id",sites[2])
    #I3=op("Id",sites[3])
    #I4=op("Id",sites[4])
    #O=Sz1*I2
    #A=Array(O,sites[1]',sites[2]',sites[1],sites[2])
    #display(A)
    #display(reshape(A,4,4))

    ←
end