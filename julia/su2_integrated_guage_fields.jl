using ITensors
using ITensorMPS
using Plots
begin
    N=20 #number of staggered fermionic sites 
        #red odd sites, green even sites 1=red_site 1 2=green site 1
    sites=siteinds("S=1/2",2*N+2);#+2 external charges
    H=OpSum()
    number_of_sites=2*N
    epsilon=1
    m=5
    r(i) = 2*i + 1
    g(i) = 2*i + 2
    g_coupling=1
    a=1
    ext1_ind=2*N+1  #index of the external charges
    ext2_ind=2*N+2
    n1=4 #postion of external charge 1 model before the jordan wigner transform
    n2=10 #" chrage 2 "
    lambda = a * g_coupling^2
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
        for j=0:N-1
            if abs(n-j)>0
                V_nm=V(n,j)
                #Q1_nQ1_N+Q2_nQ2_m
                global H+=lambda*V_nm*0.25,"S-",g(n),"S+",r(n),"S-",r(j),"S+",g(j)
                global H+=0.25*lambda*V_nm,"S-",r(n),"S+",g(n),"S-",g(j),"S+",r(j)
                #global H+=0.25*lambda*V_nm,"S-",g(n),"S+",r(n),"S-",r(j),"S+",g(j)
                #global H+=0.25*lambda*V_nm,"S-",r(n),"S+",g(n),"S-",g(j),"S+",r(j)

                #Q3_nQ3_m
                global H+=lambda/32*V_nm,"Z",g(n),"Z",g(j)
                global H+=-lambda/32*V_nm,"Z",g(n),"Z",r(j)
                global H+=-lambda/32*V_nm,"Z",r(n),"Z",g(j)
                global H+=lambda/32*V_nm,"Z",r(n),"Z",r(j)
            end
        end
    end
    for n=0:N-1
        V_n=V(n,n1)
        # interaction of vaccum with charge 1
        global H+=V_n*lambda*0.5,"Sx",ext1_ind,"S-",r(n),"Z",r(n),"S+",g(n) #X conntribution
        global H+=V_n*lambda*0.5,"Sx",ext1_ind,"S-",g(n),"Z",r(n),"S+",r(n)

        global H+=im*V_n*lambda*0.5,"Sy",ext1_ind,"S-",g(n),"S+",r(n) #Y conntribution
        global H+=-im*V_n*lambda*0.5,"Sy",ext1_ind,"S-",r(n),"S+",g(n)

        global H+=V_n*lambda*0.25,"Sz",ext1_ind,"Z",g(n)#Z conntribution
        global H+=-V_n*lambda*0.25,"Sz",ext1_ind,"Z",r(n)

            end
    
    for n=0:N-1
        V_n=V(n,n2)
        # interaction of vaccum with charge 2
        global H+=V_n*lambda*0.5,"Sx",ext2_ind,"S-",r(n),"Z",r(n),"S+",g(n) #X conntribution
        global H+=V_n*lambda*0.5,"Sx",ext2_ind,"S-",g(n),"Z",r(n),"S+",r(n)

        global H+=im*V_n*lambda*0.5,"Sy",ext2_ind,"S-",g(n),"S+",r(n) #Y conntribution
        global H+=-im*V_n*lambda*0.5,"Sy",ext2_ind,"S-",r(n),"S+",g(n)

        global H+=lambda*V_n*0.25,"Sz",ext2_ind,"Z",g(n)#Z conntribution
        global H+=-lambda*V_n*0.25,"Sz",ext2_ind,"Z",r(n)

    
    end
            #interaction between external charges
    H+=lambda*V(n1,n2)*0.5*2,"Sx",ext1_ind,"Sx",ext2_ind
    H+=lambda*V(n1,n2)*0.5*2,"Sy",ext1_ind,"Sy",ext2_ind
    H+=lambda*V(n1,n2)*0.5*2,"Sz",ext1_ind,"Sz",ext2_ind
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
