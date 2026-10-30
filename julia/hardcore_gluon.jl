using ITensors
using ITensorMPS,LinearAlgebra
using Revise
using Base.Threads
using Plots
let 
    n_matter=60
    g=10
    m=1
    gauge_link=Index(5,"GaugeLink")
    sites=Index[]
    function ITensors.siteind(
        ::SiteType"GaugeLink";
        kwargs...)
        return Index(5, "Site,GaugeLink")
    end
    println(1)
    ITensors.state(::StateName"Vaccum",::SiteType"GaugeLink")=[1.0,0,0,0]
    
    push!(sites,siteind("Electron"))
    for n=2:n_matter
        push!(sites,Index(5,"GaugeLink"))
        push!(sites,siteind("Electron"))
    end
    ITensors.op(::OpName"c_g", ::SiteType"Electron")= array(op("Cup", siteind("Electron")))
    ITensors.op(::OpName"c_r", ::SiteType"Electron")= array(op("Cdn", siteind("Electron")))

    ITensors.op(::OpName"c_g_dag", ::SiteType"Electron") = array(op("Cdagup", siteind("Electron")))
    ITensors.op(::OpName"c_r_dag", ::SiteType"Electron") = array(op("Cdagdn", siteind("Electron")))
    ITensors.op(::OpName"N_g", ::SiteType"Electron")=array(op("Nup",siteind("Electron")))
    ITensors.op(::OpName"N_r", ::SiteType"Electron")=array(op("Ndn",siteind("Electron")))
    E_sqr=Matrix(3/4*Diagonal([0 ,1, 1, 1, 1]))
    
    Urr=zeros(ComplexF64,5,5)
    Urg=zeros(ComplexF64,5,5)
    Ugr=zeros(ComplexF64,5,5)
    Ugg=zeros(ComplexF64,5,5)
    s=1/sqrt(2)
    Urr[2, 1] =  s
    Urr[1, 5] =  s
    
    Urg[3, 1] =  s
    Urg[1, 4] = -s

    Ugr[4, 1] =  s
    Ugr[1, 3] = -s

    Ugg[5, 1] =  s
    Ugg[1, 2] =  s
    ITensors.op(::OpName"Urr",::SiteType"GaugeLink") = Urr
    #@show isdefined(Main, :ITensor)
    E_sqr_op= ITensor(E_sqr,prime(gauge_link),dag(gauge_link))
    Urr_op = ITensor(Urr, prime(gauge_link), dag(gauge_link))
    Urg_op = ITensor(Urg, prime(gauge_link), dag(gauge_link))
    Ugr_op = ITensor(Ugr, prime(gauge_link), dag(gauge_link))
    Ugg_op = ITensor(Ugg, prime(gauge_link), dag(gauge_link))

    ITensors.op(::OpName"Urr", ::SiteType"GaugeLink") = Urr
    ITensors.op(::OpName"Urg", ::SiteType"GaugeLink") = Urg
    ITensors.op(::OpName"Ugr", ::SiteType"GaugeLink") = Ugr
    ITensors.op(::OpName"Ugg", ::SiteType"GaugeLink") = Ugg
    ITensors.op(::OpName"E2",  ::SiteType"GaugeLink") = E_sqr

    ITensors.op(::OpName"Urrdag", ::SiteType"GaugeLink") = adjoint(Urr)
    ITensors.op(::OpName"Urgdag", ::SiteType"GaugeLink") = adjoint(Urg)
    ITensors.op(::OpName"Ugrdag", ::SiteType"GaugeLink") = adjoint(Ugr)
    ITensors.op(::OpName"Uggdag", ::SiteType"GaugeLink") = adjoint(Ugg)
    @show hastags(sites[2], "Site")
@show hastags(sites[2], "GaugeLink")

@show size(array(op("Urr", sites[2])))
@show size(array(op("E2", sites[2])))
    H=OpSum()
    for n=1:n_matter-1
        i=2*n-1
        
        link = 2*n
        j = 2*n + 1
        H+= 1/2,"c_r_dag",i,"Urr",link,"c_r",j
        H+= 1/2,"c_r_dag",i,"Urg",link,"c_g",i+2
        H+= 1/2,"c_g_dag",i,"Ugg",link,"c_g",i+2
        H+= 1/2,"c_g_dag",i,"Ugr",link,"c_r",i+2
        
        H += 1/2, "c_r_dag", j, "Urrdag", link, "c_r", i
        H += 1/2, "c_g_dag", j, "Urgdag", link, "c_r", i
        H += 1/2, "c_g_dag", j, "Uggdag", link, "c_g", i
        H += 1/2, "c_r_dag", j, "Ugrdag", link, "c_g", i

        
        H+=g^2,"E2",i+1
        
        #H+=m*(-1)^i,"c_r_dag",i,"c_r",i,"c_g_dag",i,"c_g",i
    end

    for n in 1:n_matter
    i = 2*n - 1
    H += m*(-1)^n, "Ndn", i
    H += m*(-1)^n, "Nup", i
    end
    sweeps=40

    H_mpo=MPO(H,sites)
    psi0=random_mps(sites)
    energy,psi_0=dmrg(H_mpo,psi0,nsweeps=sweeps,maxdim=100,cutoff=1e-10)
    display(energy)
    N_r = zeros(n_matter)
    N_g = zeros(n_matter)
    for n in 1:n_matter
    i = 2*n - 1
    N_r[n] = real(expect(psi_0, "Ndn"; sites=i))
    N_g[n] = real(expect(psi_0, "Nup"; sites=i))
    end
    E=zeros(n_matter-1)
    for n in 1:n_matter-1
        i=2*n
        E[n]=real(expect(psi_0, "E2"; sites=i))
    end
        
    p1=plot(1:length(N_r),N_r,label="red",color="red")
    plot!(1:length(N_g),N_g,label="green",color="green")
    p2=plot(1:length(E),E)

     p=plot(p1, p2; layout=(2, 1))
end