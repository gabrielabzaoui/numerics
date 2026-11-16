using ITensors
using ITensorMPS
using Plots
begin
    #OBC !
    N=20 #number of fermion sites
    N_spinsites=2*N+2 # number of spin sites after jordan-wigner transform + 2 external charges
    sites=siteinds("S=1/2",N_spinsites)
    n1=4 # postion of external charges before JW transform
    n2=10
    n1_spin=2*N+1 #position of external charges on the spin chain
    n2_spin=2*N+2
    # the red (green) fermions get positiont on the odd (even) sites of the spin chain and 
    r(n)=2*n+1 #note here indexing starts at n=0
    g(n)=2*n+2
    # paramteres
    a_0=1
    m=1
    g_coupling=2
    H=OpSum()
    lambda_kin=1/(2*a_0)
    lambda_mass=m*0.5
    lambda_elec=0.5*a_0*g^2
    
    #hopping term
    for n=0:N-1
        global H+=lambda_kin,"S+",r(n+1),"Z",g(n),"S-",r(n)
        global H+=lambda_kin,"S+",g(n+1),"Z",r(n),"S-",g(n)
        
        # and the h.c.
        global H+=lambda_kin,"S-",r(n+1),"Z",g(n),"S+",r(n)
        global H+=lambda_kin,"S-",g(n+1),"Z",r(n),"S+",g(n)
    end
    
    #mass term
    for n=0:N-1
        global H+=2*lambda_mass*(-1)^n,"Id"
        global H+=-(-1)^n*lambda_mass,"Z",r(n)
        global H+=-(-1)^n*lambda_mass,"Z",g(n)
    end
    #charge operator acting on external charges
    #n1
    qx_n1=OpSum()
    qx_n1+="Sx",n1_spin
    
    qy_n1=OpSum()
    qy_n1+="Sy",n1_spin
        
    qz_n1=OpSum()
    qz_n1+="Sx",n1_spin
    
    #n2
    qx_n2=OpSum()
    qx_n2+="Sx",n2_spin
    
    qy_n2=OpSum()
    qy_n2+="Sy",n2_spin
        
    qz_n2=OpSum()
    qz_n1+="Sx",n2_spin
    
    function Q_x(n)
        op=OpSum()
        op+=-0.5,"S-",r(n),"S+",g(n)
        op+=-0.5,"S+",r(n),"S-",g(n)
        return op
    end
    function Q_y(n)
        op=OpSum()
        op+=0.5im,"S-",r(n),"S+",g(n)
        op+=-0.5im,"S+",r(n),"S-",g(n)
        return op
    end
    function Q_z(n)
        op=OpSum()
        op+=-0.25,"Z",r(n)
        op+=0.25,"Z",g(n)
        return op
    end
    function cum_Q_sum(n) #cumulative sum smaller than n
        Sx=OpSum()
        Sy=OpSum()
        Sz=OpSum()
        for i=0:n-1
            global Sx+=Q_x(i)
            global Sy+=Q_y(i)
            global Sz+=Q_z(i)
            if m==n1
                
                
                global Sx+=qx_n1

                
                global Sy+=qy_n1

                
                global Sz+=qz_n1
                
            elseif m==n2
                
                global Sx+=qx_n2

                
                global Sy+=qy_n2

                
                global Sz+=qz_n2
                
            end
        end

        return Sx,Sy,Sz
    
    end
    #electric energy 
    for n=0:N-1
        Sx,Sy,Sz=cum_Q_sum(n)
        global H+=lambda_elec*(Sx*Sx+Sy*Sy+Sz*Sz)
    
    end
    function Ln2(n)
        Lxn,Lyn,Lzn=cum_Q_sum(n)
        return Lxn*Lxn+Lyn*Lyn+Lzn*Lzn
    end
    H_mpo=MPO(H,sites)
    psi_0=random_mps(sites)
    E_0,psi=dmrg(H_mpo,psi_0;Sweeps=12;cutoff=1e-12)
    z_values = expect(psi, "Z")

    n_r = [(1-z_values[2*i+1])/2 for i in 0:N-1]
n_g = [(1-z_values[2*i+2])/2 for i in 0:N-1]
plot(n_r)
end