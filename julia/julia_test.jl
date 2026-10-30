using Plots 
using ITensors
using ITensorMPS
using Revise
using Base.Threads

let 
    N=80
    chain=siteinds("Electron",N;conserve_qns=true)
    ratio=-5:0.5:5
    
    t=1
    
    results = zeros(length(ratio))
    sweeps=12
    println("Julia threads: ", Threads.nthreads())
    println("CPU threads:   ", Sys.CPU_THREADS)
    state=[isodd(i) ? "Up" : "Dn" for i in 1:N]
    #psi=random_mps(chain,state,linkdims=10)
   
    #psi=MPS(chain,state)
   
   @time begin
    
        Threads.@threads for a in eachindex(ratio)
            
            H=OpSum()

            U=ratio[a]*t
            for i=1:N
                H+=U,"Nup",i,"Ndn",i
                if i<N
                    H+=-t,"Cdagup",i,"Cup",i+1
                    H+=-t,"Cup",i,"Cdagup",i+1
                    H+=-t,"Cdagdn",i,"Cdn",i+1
                    H+=-t,"Cdn",i,"Cdagdn",i+1
                end
                
            end
            
            H_mpo=MPO(H,chain)
            psi=MPS(chain,state)
            energy,psi0=dmrg(H_mpo,psi;nsweeps=sweeps,maxdim=100,cutoff=1e-10)
            nup_expect=expect(psi0,"Nup")
            ndn_expect=expect(psi0,"Ndn")
            #n=sum(nup_expect.+ndn_expect)/N
            double_occ = expect(psi0, "Nupdn")
            D = sum(double_occ) / N
            results[a]=D
            
        end
    end
    plot(   ratio,results,xlabel="U/t",ylabel="Average filling",marker=:circle,label="")
    #plot(filling,label="fil",ylim=(0,1))            
    
end