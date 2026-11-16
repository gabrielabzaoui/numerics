using ITensors
using ITensorMPS
using Plots
begin
    #OBC !
    N=15 #number of fermion sites
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
    g_coupling=3
    H=OpSum()
    lambda_kin=1/(2*a_0)
    lambda_mass=m*0.5
    lambda_elec=0.5*a_0*g_coupling^2
    
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
        global H+=1*lambda_mass*(-1)^n,"Id",r(n)
        global H+=1*lambda_mass*(-1)^n,"Id",g(n)
        global H+=-(-1)^n*lambda_mass,"Z",r(n)
        global H+=-(-1)^n*lambda_mass,"Z",g(n)
    end
    #charge operator acting on external charges
    #n1
    qx_n1=OpSum()
    qx_n1+=1.0,"Sx",n1_spin
    
    qy_n1=OpSum()
    qy_n1+=1.0,"Sy",n1_spin
        
    qz_n1=OpSum()
    qz_n1+=1.0,"Sz",n1_spin
    
    #n2
    qx_n2=OpSum()
    qx_n2+=1.0,"Sx",n2_spin
    
    qy_n2=OpSum()
    qy_n2+=1.0,"Sy",n2_spin
        
    qz_n2=OpSum()
    qz_n1+=1.0,"Sz",n2_spin
    
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
    function cum_Q_sum(n,seiten) #cumulative sum smaller than n
        Sx=OpSum()
        Sy=OpSum()
        Sz=OpSum()
        for i=0:n-1
             Sx+=Q_x(i)
             Sy+=Q_y(i)
             Sz+=Q_z(i)
            if i==n1
                
                
                 Sx+=qx_n1

                
                 Sy+=qy_n1

                
                 Sz+=qz_n1
                
            elseif i==n2
                
                 Sx+=qx_n2

                
                 Sy+=qy_n2

                
                 Sz+=qz_n2
                
            end
        end

        return MPO(Sx,seiten),MPO(Sy,seiten),MPO(Sz,seiten)
    
    end
    
function electric_mpo(N)
    H_E = nothing

    for n = 1:N-1
        Mx, My, Mz = cum_Q_sum(n,sites)

        Lx2 = apply(Mx, Mx; cutoff=1e-12)
        Ly2 = apply(My, My; cutoff=1e-12)
        Lz2 = apply(Mz, Mz; cutoff=1e-12)

        L2 = add(Lx2, Ly2; cutoff=1e-12)
        L2 = add(L2, Lz2; cutoff=1e-12)

        if H_E === nothing
            H_E = L2
        else
            H_E = add(H_E, L2; cutoff=1e-12) #H_E is an MPO
        end
    end

    return H_E
end
function Ln_sqr(n)
    Lx,Ly,Lz=cum_Q_sum(n,sites)
    Lx2=apply(Lx,Lx;cutoff=1e-12)
    Lz2=apply(Lz,Lz;cutoff=1e-12)
    Ly2=apply(Ly,Ly;cutoff=1e-12)
    Ln2=add(Lx2,Ly2;cutoff=1e-12)
    Ln2=add(Ln2,Lz2;cutoff=1e-12)
    return Ln2
end
   
    H_ferm=MPO(H,sites)
    H_mpo_total=add(H_ferm,electric_mpo(N),cutoff=1e-12)
    
    psi_0=random_mps(sites)
    E_0,psi=dmrg(H_mpo_total,psi_0;nsweeps=12,cutoff=1e-12,maxdim=200)
    z_values = expect(psi, "Z")

    n_r = [(1-z_values[2*i+1])/2 for i in 0:N-1]
    n_g = [(1-z_values[2*i+2])/2 for i in 0:N-1]
    L2_vals=Float64[]
    @time for i=1:N-1
        val=inner(psi,apply(Ln_sqr(i),psi))
        push!(L2_vals,real(val))
    end
    p1=scatter(n_r,color="red",label="occupation red")
    p2=scatter(n_g,color="lightgreen",label="occupation green")
    p3=plot(L2_vals,marker="circle",label="L^2")
    
    plot(p1, p2, p3, layout=(3,1), size=(700,900))
end