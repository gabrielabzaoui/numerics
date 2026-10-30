using Plots
using ITensors
using ITensorMPS
using Revise
using Base.Threads
let 
    number_of_sites=20
    sites=siteinds("S=1/2",number_of_sites)
    
    sweeps=20
    ratio=0:0.01:0.75
    J=1
    Mz_array=zeros(length(ratio))
    S_array=zeros(length(ratio))

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
    instate=["Up" for i in 1:number_of_sites]
    psi_0=MPS(sites,instate)
    @time begin 
        Threads.@threads for i in eachindex(ratio)
            
            h=ratio[i]*J
            H=OpSum()
            for i=1:number_of_sites
                H+=h,"Sx",i
                
                if i<number_of_sites
                    H+=-J,"Sz",i,"Sz",i+1
                end
            end
            if h==0
                psi_start=random_mps(sites)
            else
                psi_start=psi_0
            end
            H_mpo=MPO(H,sites)
            energy,psi=dmrg(H_mpo,psi_start;nsweeps=sweeps,maxdim=200,cutoff=1e-18)
            m=expect(psi,"Sz")
            Mz_array[i]=sum(m)/length(m)
            S_array[i]=entaglement_entropy(psi)

        
        
        





        
        end
    end
     #p1=plot(1:number_of_sites,mz,xlabel="n",ylabel="m_z",seriestype=:scatter,ylim=[-0.5,0.5],legend=false)
    p=plot(ratio,Mz_array,xlabel="h/J",ylabel="Mz",label="Mz" ,legend=:topleft,ylim=[-0.1,0.5])
    plot!(
    twinx(p),
    ratio, S_array,
    ylabel = "S",
    label = "S",
    color = :red,
    legend = :topright
)
    
    #display(plot(p1, p2; layout=(2, 1)))
    savefig(p,"phasetransition.png")
    display(p)

    
end