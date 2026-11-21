using ITensors,ITensorMPS,HDF5,Plots
begin
    files = readdir("res/")
    filename="res/su2_res_4_16.h"
    i=1
    for filename in files
    
    h5open("res/$(filename)","r") do f
    #h5open(filename,"r") do f
    println(i)
       # for key in keys(f)
        #    println(key)
        #end
        psi = read(f, "psi",MPS)
     #   println(typeof(psi))
    #println(inner(psi',psi))
    n1=read(f,"n1")
    n2=read(f,"n2")
    Qx2 = read(f, "Qx2",MPO)
    Qy2 = read(f, "Qy2",MPO)
    Qz2 = read(f, "Qz2",MPO)

    qx2 = real(inner(psi', Qx2, psi))
    qy2 = real(inner(psi', Qy2, psi))
    qz2 = real(inner(psi', Qz2, psi))

    Qtotal2 = qx2 + qy2 + qz2
    L2_vals=read(f,"L2_vals")
    n_r=read(f,"n_r")
    n_g=read(f,"n_g")
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
    
    savefig(p,"res_plot/$(splitext(filename)[1]).pdf")
    println("<Qx²> = ", qx2)
    println("<Qy²> = ", qy2)
    println("<Qz²> = ", qz2)
    println("<Qtotal²> = ", Qtotal2)
    println("n_r")
    println(n_r)
    println("n_g")
    println(n_g)
    println("L^2")
    println(L2_vals)    
    end
    global i+=1
end
end
