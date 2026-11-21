using ITensors
using ITensorMPS
using CUDA,HDF5
#using Plots
begin
    N=20 #number of staggered fermionic sites 
        #red odd sites, green even sites 1=red_site 1 2=green site 1
    sites=siteinds("S=1/2",2*N+2);#+2 external charges
    H=OpSum()
    number_of_sites=2*N
    
    m=1
    r(i) = 2*i + 1
    g(i) = 2*i + 2
    g_coupling=2
    a=1
    epsilon=1/(2*a)
    ext1_ind=2*N+1  #index of the external charges
    ext2_ind=2*N+2
    n1=parse(Int,ARGS[1])#postion of external charge 1 model before the jordan wigner transform
    n2=16 #" chrage 2 "
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
function total_colour_operators(Itensor_sites)

    Qx_os = OpSum()
    Qy_os = OpSum()
    Qz_os = OpSum()

    for i = 0:N-1

        Qx_os +=  0.5, "S-", r(i), "S+", g(i)
        Qx_os +=  0.5, "S+", r(i), "S-", g(i)

        Qy_os += -0.5im, "S-", r(i), "S+", g(i)
        Qy_os +=  0.5im, "S+", r(i), "S-", g(i)

        Qz_os +=  0.25, "Z", g(i)
        Qz_os += -0.25, "Z", r(i)

    end

    Qx_os += 1.0, "Sx", ext1_ind
    Qx_os += 1.0, "Sx", ext2_ind

    Qy_os += 1.0, "Sy", ext1_ind
    Qy_os += 1.0, "Sy", ext2_ind

    Qz_os += 1.0, "Sz", ext1_ind
    Qz_os += 1.0, "Sz", ext2_ind

    return MPO(Qx_os, Itensor_sites),
           MPO(Qy_os, Itensor_sites),
           MPO(Qz_os, Itensor_sites)
end

# ============================================
# Construct total colour operators
# ============================================

Qx, Qy, Qz = total_colour_operators(sites)

# ============================================
# Construct total Casimir
# ============================================

cutoff_penalty = 1e-14

Qx2_mpo = apply(Qx, Qx;
                cutoff=cutoff_penalty,
                maxdim=4000)

Qy2_mpo = apply(Qy, Qy;
                cutoff=cutoff_penalty,
                maxdim=4000)

Qz2_mpo = apply(Qz, Qz;
                cutoff=cutoff_penalty,
                maxdim=4000)

Qtotal2_mpo = add(Qx2_mpo, Qy2_mpo;
                  cutoff=cutoff_penalty,
                  maxdim=4000)

Qtotal2_mpo = add(Qtotal2_mpo, Qz2_mpo;
                  cutoff=cutoff_penalty,
                  maxdim=4000)

# ============================================
# Add colour neutrality penalty
# ============================================

kappa = 10.0

H_base = MPO(H, sites)

H_penalty = kappa * Qtotal2_mpo

H_mpo = add(H_base, H_penalty;
            cutoff=cutoff_penalty,
            maxdim=4000)

# ============================================
# DMRG
# ============================================

psi_0 = random_mps(sites)
H_gpu=cu(H_mpo)
psi_0gpu=cu(psi_0)
E_0, psi_gpu = dmrg(H_gpu, psi_0gpu;
                nsweeps=100,
                cutoff=1e-14,
                maxdim=800)
CUDA.synchronize()
psi=cpu(psi_gpu)
# ============================================
# Measure total colour
# ============================================

Qx2 = real(inner(Qx, psi, Qx, psi))
Qy2 = real(inner(Qy, psi, Qy, psi))
Qz2 = real(inner(Qz, psi, Qz, psi))

Qtotal2 = Qx2 + Qy2 + Qz2

println("<Qx²> = ", Qx2)
println("<Qy²> = ", Qy2)
println("<Qz²> = ", Qz2)
println("<Qtotal²> = ", Qtotal2)

# ============================================
# Matter occupation
# ============================================

z_values = expect(psi, "Z")

n_r = [(1-z_values[2*i+1])/2 for i in 0:N-1]
n_g = [(1-z_values[2*i+2])/2 for i in 0:N-1]
"""
#scatter(1:N, n_r,
        marker=:square,
        color="red",
        label="n_r")

#scatter!(1:N, n_g,
         marker=:circle,
         color="lightgreen",
         label="n_g")
"""

# ============================================
# Electric flux profile
# ============================================

function cumulative_colour_operators(sites, x)

    Qx_os = OpSum()
    Qy_os = OpSum()
    Qz_os = OpSum()

    # Dynamical matter charges to the left of the cut
    for i in 0:x

        Qx_os +=  0.5, "S-", r(i), "S+", g(i)
        Qx_os +=  0.5, "S+", r(i), "S-", g(i)

        Qy_os += -0.5im, "S-", r(i), "S+", g(i)
        Qy_os +=  0.5im, "S+", r(i), "S-", g(i)

        Qz_os +=  0.25, "Z", g(i)
        Qz_os += -0.25, "Z", r(i)

    end

    # Include external charges once they are to the left of the cut
    if n1 <= x
        Qx_os += 1.0, "Sx", ext1_ind
        Qy_os += 1.0, "Sy", ext1_ind
        Qz_os += 1.0, "Sz", ext1_ind
    end

    if n2 <= x
        Qx_os += 1.0, "Sx", ext2_ind
        Qy_os += 1.0, "Sy", ext2_ind
        Qz_os += 1.0, "Sz", ext2_ind
    end

    return MPO(Qx_os, sites),
           MPO(Qy_os, sites),
           MPO(Qz_os, sites)
end


# Calculate <L_x^2 + L_y^2 + L_z^2> at each cut

flux_squared = zeros(Float64, N)
Lx_op=MPO[]
Ly_op=MPO[]
Lz_op=MPO[]
for x in 0:N-1

    Lx, Ly, Lz = cumulative_colour_operators(sites, x)

    Fx = real(inner(Lx, psi, Lx, psi))
    Fy = real(inner(Ly, psi, Ly, psi))
    Fz = real(inner(Lz, psi, Lz, psi))

    flux_squared[x+1] = Fx + Fy + Fz
    push!(Lx_op,Lx)
    push!(Ly_op,Ly)
    push!(Lz_op,Lz)
end
output_file = "SU2_N$(N)_n1$(n1)_n2$(n2).h5"

h5open(output_file, "w") do f

    # Ground state
    write(f, "psi", psi)

    # Energy
    write(f, "E0", E_0)

    # Matter occupation
    write(f, "n_r", n_r)
    write(f, "n_g", n_g)

    # Electric flux
    write(f, "flux_squared", flux_squared)

    # Total colour charge
    write(f, "Qx2", Qx2)
    write(f, "Qy2", Qy2)
    write(f, "Qz2", Qz2)
    write(f, "Qtotal2", Qtotal2)

    # Physical parameters
    write(f, "N", N)
    write(f, "n1", n1)
    write(f, "n2", n2)
    
    write(f, "m", m)
    write(f, "g_coupling", g_coupling)
    write(f, "a", a)
    write(f, "lambda", lambda)
    write(f, "kappa", kappa)
    for x in 1:N
        write(f, "Lx_$x", Lx_op[x])
        write(f, "Ly_$x", Ly_op[x])
        write(f, "Lz_$x", Lz_op[x])
    end
    # DMRG parameters
    write(f, "nsweeps", 100)
write(f, "maxdim", 800)
write(f, "cutoff", 1e-14)
    """
    write(f, "nsweeps", 12)
    write(f, "maxdim", 200)
    write(f, "cutoff", 1e-12)
    """
    # Mapping information
    write(f, "MPS_ordering", "r(0),g(0),r(1),g(1),...,external1,external2")
end

println("Saved results to ", output_file)
# Plot the electric flux
#=
p1=scatter(n_r,color="red",label="occupation red")
    p2=scatter(n_g,color="lightgreen",label="occupation green")
    p3=plot(flux_squared,marker="circle",label="L^2")
    
    
# Mark external charge positions
vline!([n1, n2],
       linestyle=:dash,
       linewidth=2,
       label="External charges")
plot(p1, p2, p3, layout=(3,1), size=(700,900))
=#
end
