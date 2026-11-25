using ITensors, ITensorMPS
using Plots

red(n)=2*n+1
green(n)=2*n+2
function get_kinetic_hamiltonian(epsilon,N,seiten)
    os=OpSum()
    #number_of_sites=2*N
    for i=0:N-2
        os+=epsilon,"S+",red(i),"Z",green(i),"S-",red(i+1)
        os+=epsilon,"S-",red(i),"Z",green(i),"S+",red(i+1) #hermitian conjugate

        os+=epsilon,"S+",green(i),"Z",red(i+1),"S-",green(i+1)
        os+=epsilon,"S-",green(i),"Z",red(i+1),"S+",green(i+1) #hermitian conjugate


    end
    return os
end

function get_mass_term(m,N,seiten)
    #number_of_sites=2*N
    os=OpSum()
    for i=0:N-1
        os+=m*(-1)^i,"Z",red(i)
        os+=m*(-1)^i,"Z",green(i)
        os+=2*m*(-1)^i,"Id",red(i)
        os+=2*m*(-1)^i,"Id",green(i)

    end
    return os 
end

Qx(n) = let os = OpSum()
    os += (-1im/2, "S+", red(n), "S-", green(n))
    os += ( 1im/2, "S-", red(n), "S+", green(n))
    os 
end

Qy(n) = let os = OpSum()
    os += (-1/2, "S+", red(n), "S-", green(n))
    os += ( -1/2, "S-", red(n), "S+", green(n))
    os
end
Qz(n) = let os = OpSum()
    os += (1/4, "Z", red(n))
    os += ( -1/4, "Z", green(n))
    os 
end
V(n,m)=-0.5*abs(n-m)
function get_electric_hamiltonian(g,a,N,n1,n2)
    os=OpSum()
    
    for n=0:N-1
        for m=0:N-1
            
            if n==m 
                continue

            else 
                pref = 0.5 * g^2 * a * V(n,m)
                #QxnQxm+QynQym
                os += pref/2, "S+", red(n), "S-", green(n),
                 "S-", red(m), "S+", green(m)

                os += pref/2, "S-", red(n), "S+", green(n),
                 "S+", red(m), "S-", green(m)

                #QznQzm
                os +=  pref/16, "Z", red(n), "Z", red(m)
                os += -pref/16, "Z", red(n), "Z", green(m)
                os += -pref/16, "Z", green(n), "Z", red(m)
                os +=  pref/16, "Z", green(n), "Z", green(m)
            end
        end
        #interacting with external charges
        Vnn1=V(n,n1)
        Vnn2=V(n,n2)
        factor=g^2*a
        #x
        os+=factor*Vnn1*(-im/2), "S+", red(n), "S-", green(n),"Sx",2*N+1
        os+=factor*Vnn1*(im/2), "S-", red(n), "S+", green(n),"Sx",2*N+1

        os+=factor*Vnn2*(-im/2), "S+", red(n), "S-", green(n),"Sx",2*N+2
        os+=factor*Vnn2*(im/2), "S-", red(n), "S+", green(n),"Sx",2*N+2
        #y
        os+=factor*Vnn1*(-1/2),"S+", red(n), "S-", green(n),"Sy",2*N+1
        os+=factor*Vnn1*(-1/2),"S-", red(n), "S+", green(n),"Sy",2*N+1
        
        os+=factor*Vnn2*(-1/2),"S+", red(n), "S-", green(n),"Sy",2*N+2
        os+=factor*Vnn2*(-1/2),"S-", red(n), "S+", green(n),"Sy",2*N+2
        #z 
        os+=factor*Vnn1*(1/4),"Z",red(n),"Sz",2*N+1
        os+=factor*Vnn1*(-1/4),"Z",green(n),"Sz",2*N+1

        os+=factor*Vnn2*(1/4),"Z",red(n),"Sz",2*N+2
        os+=factor*Vnn2*(-1/4),"Z",green(n),"Sz",2*N+2
    end
    Vn1n2=V(n1,n2)

    os+=Vn1n2*g^2*a,"Sx",2*N+1,"Sx",2*N+2
    os+=Vn1n2*g^2*a,"Sy",2*N+1,"Sy",2*N+2
    os+=Vn1n2*g^2*a,"Sz",2*N+1,"Sz",2*N+2
    
    return os
end
function get_penalty_mpo(lambda,N,seiten)
    Qx_total = OpSum()
    Qy_total = OpSum()
    Qz_total = OpSum()

    for n in 0:N-1
        Qx_total += -im/2, "S+", red(n), "S-", green(n)
        Qx_total +=  im/2, "S-", red(n), "S+", green(n)

        Qy_total += -1/2, "S+", red(n), "S-", green(n)
        Qy_total += -1/2, "S-", red(n), "S+", green(n)

        Qz_total += 1/4, "Z", red(n)
        Qz_total += -1/4, "Z", green(n)
    end

    # Add external charges
    for j in (2*N+1, 2*N+2)
        Qx_total += 1.0, "Sx", j
        Qy_total += 1.0, "Sy", j
        Qz_total += 1.0, "Sz", j
    end
    Qx_total_mpo=MPO(Qx_total,seiten)
    Qy_total_mpo=MPO(Qy_total,seiten)
    Qz_total_mpo=MPO(Qz_total,seiten)
    Mx=apply(Qx_total_mpo,Qx_total_mpo; cutoff=1e-12)
    My=apply(Qy_total_mpo,Qy_total_mpo; cutoff=1e-12)
    Mz=apply(Qz_total_mpo,Qz_total_mpo; cutoff=1e-12)
    H_penalty=add(Mx,My)
    H_penalty=lambda*add(H_penalty,Mz)
    return H_penalty
end
function get_hamiltonian_mpo(epsilon,m,g,a,N,n1,n2,lambda,seiten)
    kin=get_kinetic_hamiltonian(epsilon,N,seiten)
    mass=get_mass_term(m,N,seiten)
    elec=get_electric_hamiltonian(g,a,N,n1,n2)
    H_penalty=get_penalty_mpo(lambda,N,seiten)
    H_opsum=kin+mass+elec
    H_total_mpo=add(MPO(H_opsum,seiten),H_penalty)

    return H_total_mpo
    
end
function get_flux_mpo(n, N, sites;n1=-1, n2=-1)

    Ex = OpSum()
    Ey = OpSum()
    Ez = OpSum()

    # Dynamical charges to the left of and including n
    for m in 0:n

        Ex += -im/2, "S+", red(m), "S-", green(m)
        Ex +=  im/2, "S-", red(m), "S+", green(m)

        Ey += -1/2, "S+", red(m), "S-", green(m)
        Ey += -1/2, "S-", red(m), "S+", green(m)

        Ez += 1/4, "Z", red(m)
        Ez += -1/4, "Z", green(m)
    end

    # Include external charges if they are to the left of the link
    if n1 >= 0 && n1 <= n
        Ex += 1.0, "Sx", 2*N+1
        Ey += 1.0, "Sy", 2*N+1
        Ez += 1.0, "Sz", 2*N+1
    end

    if n2 >= 0 && n2 <= n
        Ex += 1.0, "Sx", 2*N+2
        Ey += 1.0, "Sy", 2*N+2
        Ez += 1.0, "Sz", 2*N+2
    end

    return MPO(Ex, sites), MPO(Ey, sites), MPO(Ez, sites)
end
function get_electric_flux(psi,N,sites;n1=-1,n2=-1)

    flux = zeros(N)

    for n in 0:N-1

        Ex, Ey, Ez = get_flux_mpo(n,N,sites;n1=n1,n2=n2)
        #=
        Ex2 = apply(Ex, Ex; cutoff=1e-12)
        Ey2 = apply(Ey, Ey; cutoff=1e-12)
        Ez2 = apply(Ez, Ez; cutoff=1e-12)=#

        flux[n+1] = real(
            inner(Ex,psi, Ex, psi) +
            inner(Ey,psi,Ey, psi) +
            inner(Ez,psi, Ez, psi)
        )
    end

    return flux
end
function Start_simulation(epsilon,m,g,a,N,n1,n2,lambda;calc::Bool=true,plt::Bool=true,nsweeps=15,cutoff=1e-12,maxdim=400)
    sites=siteinds("S=1/2",2*N+2)
    H_mpo=get_hamiltonian_mpo(epsilon,m,g,a,N,n1,n2,lambda,sites)
    psi_0=random_mps(sites)
    energy,psi=dmrg(H_mpo,psi_0;cutoff=cutoff,nsweeps=nsweeps,maxdim=maxdim)
    if calc
        
        Z=expect(psi,"Z")
        n_r = [(1+Z[2*i+1])/2 for i in 0:N-1]
        n_g = [(1+Z[2*i+2])/2 for i in 0:N-1]
        flux=get_electric_flux(psi,N,sites;n1=n1,n2=n2)
        if plt
            p1=scatter(0:N-1, n_r, label="Red",color="red")
            ylims!(p1,0,1)
            p2=scatter(0:N-1, n_g, label="Green",color="lightgreen")
            ylims!(p2,0,1)
            p3=plot( flux,
                label="Electric flux",
                xlabel="Link",
                ylabel="⟨E²⟩",marker=:circle)
                ylims!(p3,0,1.0)
            p_tot=plot(p1, p2,p3, layout=(3,1), size=(700,900))
            display(p_tot)
        end
        return energy, psi,n_r,n_g, flux

    end
    return energy, psi
end