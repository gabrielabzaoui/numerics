using ITensors
using ITensorMPS
using Revise
using Base.Threads
using Plots
let 
    number_of_sites=30
    sites=siteinds("S=1/2",number_of_sites)
    H=OpSum()
    sweeps=30
    dt=0.02
    time_steps=40
    Mz_array=zeros(time_steps)
    S_arrray=zeros(time_steps)

    h=20
    J=-8
    for i=1:number_of_sites
        H+=h,"Sx",i
        
        if i<number_of_sites
            H+=-J,"Sz",i,"Sz",i+1
        end
        
    
    
    
    end
    function entaglement_entropy(psi_mps)
        b=length(psi_mps)÷2
        psi = orthogonalize(psi_mps, b)
        U,S,V = svd(psi[b], (linkinds(psi, b-1)..., siteinds(psi, b)...))
        SvN = 0.0
        for n=1:dim(S, 1)
            p = S[n,n]^2
            SvN -= p * log(p)
        end
        return SvN        
    end
    H_mpo=MPO(H,sites)
    instate=["Up" for i in 1:number_of_sites]
    #psi0=random_mps(sites)
    psi0=MPS(sites,instate)
    energy,psi=dmrg(H_mpo,psi0,nsweeps=sweeps,maxdim=100,cutoff=1e-10)
    mz=expect(psi,"Sz")
    Mz=sum(mz)/length(mz)
    println(linkdims(psi0))
    @time for i=1:time_steps
        psi0 = tdvp(H_mpo,-1im*dt,psi0;nsteps=1, time_step=-1im * dt,normalize=false,cutoff=1e-13,nsite=1)
        m=expect(psi0,"Sz")
        Mz_array[i]=sum(m)/length(m)
        S_arrray[i]=entaglement_entropy(psi0)
        #psi0=psi2
        println(linkdims(psi0))


    end
    times=1:time_steps
    times=dt*times
    p1=plot(1:number_of_sites,mz,xlabel="n",ylabel="m_z",seriestype=:scatter,ylim=[-0.5,0.5],legend=false)
    p2=plot(times,[Mz_array ],xlabel="t",ylabel="Mz",label=["Mz" ],legend=:topleft,ylim=[-0.5,0.5])
    plot!(
    twinx(p2),
    times, S_arrray,
    ylabel = "S",
    label = "S",
    color = :red,
    legend = :topright
)
    
    
    p=plot(p1, p2; layout=(2, 1))
    #savefig(p,"entanglement.png")
    #println(Mz)
end