using ITensors
using ITensorMPS
using Plots
begin
    N=20 #number of staggered fermionic sites 
        #red odd sites, green even sites 1=red_site 1 2=green site 1
    sites=siteinds("S=1/2",2*N);
    H=OpSum()
    number_of_sites=2*N
    epsilon=1
    m=5
    r(i) = 2*i + 1
    g(i) = 2*i + 2
    g_coupling=1
    a=1
    function number_r(sites, i)
        os = OpSum()
        os += 0.5, "Id", r(i)
        os += -0.5, "Z", r(i)
        return MPO(os, sites)
    end

    function number_g(sites, i)
        os = OpSum()
        os += 0.5, "Id", g(i)
        os += -0.5, "Z", g(i)
        return MPO(os, sites)
    end
    function V(n,m)
        return -0.5*abs(n-m)
    end
    for i=0:N-2
        r_i   = 2*i + 1
        g_i   = 2*i + 2
        r_ip1 = 2*(i+1) + 1
        g_ip1 = 2*(i+1) + 2
        global H+=epsilon,"S-",r_i,"Z",g_i,"S+",r_ip1 #ψ^†_rψ_r
        global H+=epsilon,"S+",r_i,"Z",g_i,"S-",r_ip1 #h.c

        global H+=epsilon,"S-",g_i,"Z",r_ip1,"S+",g_ip1 #ψ^†_gψ_g
        global H+=epsilon,"S+",g_i,"Z",r_ip1,"S-",g_ip1 #h.c

        


    end
    for i=0:N-1
        r_i   = 2*i + 1
        g_i   = 2*i + 2
        r_ip1 = 2*(i+1) + 1
        g_ip1 = 2*(i+1) + 2
        #mass term
        global H+=(-1)^(i+1)*m*0.5,"Z",g_i
        global H+=(-1)^(i+1)*m*0.5,"Z",r_i
        global H+=(-1)^(i)*m*0.5,"Id",g_i
        global H+=(-1)^(i)*m*0.5,"Id",r_i
    end
    #H_electric
    for n=0:N-1
        for m=0:N-1
            if abs(n-m)>0
                #Q1_nQ1_M+Q2_nQ2_m
                global H+=V(n,m)*0.5,"S-",g(n),"S+",r(n),"S-",r(m),"S+",g(m)
                global H+=V(n,m)*0.5,"S-",r(n),"S+",g(n),"S-",g(m),"S+",r(m)

                #Q3_nQ3_m
                global H+=0.5*a*(g_coupling^2)*V(n,m)*0.0625,"Z",g(n),"Z",g(m)
                global H+=0.5*a*(g_coupling^2)*V(n,m)*-0.0625,"Z",g(n),"Z",r(m)
                global H+=0.5*a*(g_coupling^2)*V(n,m)*-0.0625,"Z",r(n),"Z",g(m)
                global H+=0.5*a*(g_coupling^2)*V(n,m)*0.0625,"Z",r(n),"Z",r(m)
            end
        end
    end
    d=0.7
    H_mpo=MPO(H,sites);
    psi_0=random_mps(sites)
    E_0,psi=dmrg(H_mpo,psi_0;nsweeps=10,cutoff=1e-12,maxdim=200)
    z_values = expect(psi, "Z")
    n_r = [(1 - z_values[2*i+1])/2 for i in 0:N-1]
    n_g = [(1 - z_values[2*i+2])/2 for i in 0:N-1]
    scatter(1:length(n_r),n_r,marker=:square,color="red",label="n_r")
    scatter!(1:length(n_g),n_g,marker=:circle,color="lightgreen",label="n_g")
end
