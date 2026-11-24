using ITensors,ITensorMPS,HDF5,Plots
begin

    filename="/home/gabriel/Universität/LMU/sim_res/SU2_N20_n115_n216.h5"
    files = readdir("/home/gabriel/Universität/LMU/sim_res")
    for filename in files
        h5open("/home/gabriel/Universität/LMU/sim_res/$filename","r") do f 
            
            n_r=read(f,"n_r")
            n_g=read(f,"n_g")
            L2_vals=read(f,"flux_squared")
            n1=read(f,"n1")
            n2=read(f,"n2")
            psi=read(f,"psi",MPS)
            
            Qtotal2=read(f,"Qtotal2")
            #Qtotal2=inner(psi',Qtotal2,psi)
            
            p1=scatter(n_r,label="n_r",color=:red,marker=:circle)
            p2=scatter(n_g,label="n_g",color=:lightgreen,marker=:circle)

            p3=plot(L2_vals,marker=:circle,label="L^2")
            vline!([n1, n2],
            linestyle=:dash,
            linewidth=2,
            label="External charges")
            
            annotate!(2, maximum(L2_vals),
            text("⟨Q²_total⟩ = $(round(Qtotal2, digits=6))", 10, :left))
            
            p=plot(p1, p2, p3, layout=(3,1), size=(700,900))
            println("L2=")
            println(L2_vals)
            println("n_r=")
            println(n_r)
            
            #savefig(p,"res2_plot/n1=$n1.pdf")

            
        end
    end 
   
    
    
end