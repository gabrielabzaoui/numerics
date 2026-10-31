using ITensors
using ITensorMPS
using Revise
using Base.Threads
using Plots
using LinearAlgebra
let

    number_of_sites=30
    sites=siteinds("S=1/2",number_of_sites)
    H=OpSum()
    sweeps=30
    dt=0.01
    time_steps=50
    Mz_array_pad=zeros(time_steps)
    S_arrray_pad=zeros(time_steps)
    Mz_array=zeros(time_steps)
    S_arrray=zeros(time_steps)
    norm=zeros(time_steps)
    norm_pad=zeros(time_steps)
    h=20
    J=-8
    for i=1:number_of_sites
        H+=h,"Sx",i
        
        if i<number_of_sites
            H+=-J,"Sz",i,"Sz",i+1
        end
        
    
    
    
    end
    H_mpo=MPO(H,sites)
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
    #psi0=random_mps(sites)
    ψ_0 =MPS(sites,instate)
    ψ_pad= 0.0 * random_mps(sites; linkdims=50)

    ψ_pad = +(ψ_0, ψ_pad; alg="directsum")
    @time for i=1:time_steps
        
       
        task_one= Threads.@spawn begin
            ψ_0 = tdvp(H_mpo,-1im*dt,ψ_0;nsteps=1, time_step=-1im * dt,normalize=false,cutoff=1e-13,nsite=2)
            m = expect(ψ_0, "Sz")
            Mz_array[i] = sum(m) / length(m)
            S_arrray[i] = entaglement_entropy(ψ_0)
            norm[i]=sqrt(real(inner(ψ_0,ψ_0)))

        end
        
    
        task_two = Threads.@spawn begin
            
        
            ψ_pad = tdvp(H_mpo,-1im*dt,ψ_pad;nsteps=1, time_step=-1im * dt,normalize=false,cutoff=1e-13,nsite=1)
            m = expect(ψ_pad, "Sz")
            Mz_array_pad[i] = sum(m) / length(m)
            S_arrray_pad[i] = entaglement_entropy(ψ_pad)
            norm_pad[i]=sqrt(real(inner(ψ_pad,ψ_pad)))
            
        end
        #ψ_0, Mz_array[i], S_arrray[i] = fetch(task_one)
        #ψ_pad, Mz_array_pad[i], S_arrray_pad[i] = fetch(task_two)
        
        wait(task_one)
        wait(task_two)
        println(i)
        @show linkdims(ψ_0)
        #println( "Minimale norm:", minimum(norm_pad))
    end

    times=1:time_steps
    times=dt*times
    p1=plot(times,[Mz_array ],xlabel="t",ylabel="Mz",label="Mz" ,legend=:topleft,ylim=[-0.5,0.5])
    plot!(
    twinx(p1),
    times, S_arrray,
    ylabel = "S",
    label = "S",
    color = :red,
    legend = :topright)
    p2=plot(times,[Mz_array_pad ],xlabel="t",ylabel="Mz",label="Mz" ,legend=:topleft,ylim=[-0.5,0.5])
    plot!(
    twinx(p2),
    times, S_arrray_pad,
    ylabel = "S",
    label = "S",
    color = :red,
    legend = :topright)
    
    p3=plot(times,[norm norm_pad],xlabel="t",ylabel="",label=[ "no padding" "with padding"],ylims=(minimum(norm_pad),1))
    p=plot(p1, p2,p3; layout=(3, 1))
    
    savefig(p,"padding.pdf")
end