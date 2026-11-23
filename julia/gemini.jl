using ITensors
using ITensorMPS
using HDF5
using Plots

function run_simulation(args=ARGS)
    # ============================================
    # Parameters
    # ============================================
    N = 20                  # Number of staggered fermionic sites (0 to N-1)
    m = 1.0                 # Fermion mass
    g_coupling = 2.0        # Gauge coupling
    a = 1.0                 # Lattice spacing
    epsilon = 1.0 / (2 * a)
    lambda = a * g_coupling^2
    kappa = 20.0            # Penalty weight to enforce total SU(2) singlet

    # External static charge positions (0-based staggered site indices)
    n1 = isempty(args) ? 4 : parse(Int, args[1])
    n2 = 16
    @assert 0 <= n1 < N && 0 <= n2 < N "External charges must lie in 0:(N-1)"

    # Site indexing: 1..2N are dynamical matter (r(i), g(i)), 2N+1..2N+2 are static charges
    r(i) = 2 * i + 1
    g(i) = 2 * i + 2
    ext1_ind = 2 * N + 1
    ext2_ind = 2 * N + 2

    sites = siteinds("S=1/2", 2 * N + 2)

    # 1D Coulomb kernel for open boundary conditions
    V(n, j) = -0.5 * abs(n - j)

    # Helper to add local color charge operators Q^a_n to an OpSum
    function add_matter_charge!(os::OpSum, coeff::Number, comp::Symbol, n::Int)
        if comp == :x
            os += 0.5 * coeff, "S-", r(n), "S+", g(n)
            os += 0.5 * coeff, "S+", r(n), "S-", g(n)
        elseif comp == :y
            os += -0.5im * coeff, "S-", r(n), "S+", g(n)
            os +=  0.5im * coeff, "S+", r(n), "S-", g(n)
        elseif comp == :z
            os +=  0.25 * coeff, "Z", g(n)
            os += -0.25 * coeff, "Z", r(n)
        end
    end

    # ============================================
    # Build Hamiltonian OpSum
    # ============================================
    H = OpSum()

    # 1. Staggered kinetic (hopping) term
    for i in 0:(N - 2)
        r_i, g_i     = r(i), g(i)
        r_ip1, g_ip1 = r(i + 1), g(i + 1)

        H += epsilon, "S-", r_i, "Z", g_i, "S+", r_ip1
        H += epsilon, "S+", r_i, "Z", g_i, "S-", r_ip1

        H += epsilon, "S-", g_i, "Z", r_ip1, "S+", g_ip1
        H += epsilon, "S+", g_i, "Z", r_ip1, "S-", g_ip1
    end

    # 2. Staggered mass term: m * (-1)^i * (n_r(i) + n_g(i)), where n = (I - Z)/2
    for i in 0:(N - 1)
        sign_i = (-1)^i
        H += -0.5 * m * sign_i, "Z", r(i)
        H += -0.5 * m * sign_i, "Z", g(i)
        H +=  0.5 * m * sign_i, "Id", r(i)
        H +=  0.5 * m * sign_i, "Id", g(i)
    end

    # 3. Electric energy: H_E = (lambda / 2) * sum_{n=0}^{N-2} L_n^2
    for n in 0:(N - 1)
        w_self = 0.5 * lambda * (N - 1 - n)
        if w_self != 0
            H +=  0.375 * w_self, "Id", r(n)
            H += -0.375 * w_self, "Z", r(n), "Z", g(n)
        end

        for j in 0:(N - 1)
            if n != j
                c_nj = 0.5 * lambda * V(n, j)
                H += 0.5 * c_nj, "S-", r(n), "S+", g(n), "S+", r(j), "S-", g(j)
                H += 0.5 * c_nj, "S+", r(n), "S-", g(n), "S-", r(j), "S+", g(j)

                H +=  (1.0 / 16.0) * c_nj, "Z", g(n), "Z", g(j)
                H += -(1.0 / 16.0) * c_nj, "Z", g(n), "Z", r(j)
                H += -(1.0 / 16.0) * c_nj, "Z", r(n), "Z", g(j)
                H +=  (1.0 / 16.0) * c_nj, "Z", r(n), "Z", r(j)
            end
        end

        for (ext_pos, ext_idx) in ((n1, ext1_ind), (n2, ext2_ind))
            c_ext = lambda * V(n, ext_pos)
            if c_ext != 0
                H +=  0.5 * c_ext,   "Sx", ext_idx, "S-", r(n), "S+", g(n)
                H +=  0.5 * c_ext,   "Sx", ext_idx, "S+", r(n), "S-", g(n)

                H += -0.5im * c_ext, "Sy", ext_idx, "S-", r(n), "S+", g(n)
                H +=  0.5im * c_ext, "Sy", ext_idx, "S+", r(n), "S-", g(n)

                H +=  0.25 * c_ext,  "Sz", ext_idx, "Z", g(n)
                H += -0.25 * c_ext,  "Sz", ext_idx, "Z", r(n)
            end
        end
    end

    # 4. External charge - external charge interaction
    c_12 = lambda * V(n1, n2)
    H += c_12, "Sx", ext1_ind, "Sx", ext2_ind
    H += c_12, "Sy", ext1_ind, "Sy", ext2_ind
    H += c_12, "Sz", ext1_ind, "Sz", ext2_ind

    # ============================================
    # Total Color Operators & Singlet Penalty
    # ============================================
    function cumulative_colour_operators(sites, x_max::Int)
        Qx_os, Qy_os, Qz_os = OpSum(), OpSum(), OpSum()

        for i in 0:x_max
            add_matter_charge!(Qx_os, 1.0, :x, i)
            add_matter_charge!(Qy_os, 1.0, :y, i)
            add_matter_charge!(Qz_os, 1.0, :z, i)
        end

        if n1 <= x_max
            Qx_os += 1.0, "Sx", ext1_ind
            Qy_os += 1.0, "Sy", ext1_ind
            Qz_os += 1.0, "Sz", ext1_ind
        end

        if n2 <= x_max
            Qx_os += 1.0, "Sx", ext2_ind
            Qy_os += 1.0, "Sy", ext2_ind
            Qz_os += 1.0, "Sz", ext2_ind
        end

        return MPO(Qx_os, sites), MPO(Qy_os, sites), MPO(Qz_os, sites)
    end

    Qx, Qy, Qz = cumulative_colour_operators(sites, N - 1)

    cutoff_mpo = 1e-12
    Qx2_mpo = apply(Qx, Qx; cutoff=cutoff_mpo, maxdim=500)
    Qy2_mpo = apply(Qy, Qy; cutoff=cutoff_mpo, maxdim=500)
    Qz2_mpo = apply(Qz, Qz; cutoff=cutoff_mpo, maxdim=500)

    Qtotal2_mpo = add(Qx2_mpo, Qy2_mpo; cutoff=cutoff_mpo, maxdim=500)
    Qtotal2_mpo = add(Qtotal2_mpo, Qz2_mpo; cutoff=cutoff_mpo, maxdim=500)

    H_base = MPO(H, sites)
    H_mpo = add(H_base, kappa * Qtotal2_mpo; cutoff=cutoff_mpo, maxdim=1000)

    # ============================================
    # Initial State (Strong-Coupling Vacuum + Noise) & CPU DMRG
    # ============================================
    init_states = [isodd(div(k - 1, 2)) ? "Dn" : "Up" for k in 1:(2 * N)]
    append!(init_states, ["Up", "Dn"])
    psi_0 = random_mps(sites, init_states; linkdims=10)

    sweeps = Sweeps(30)
    setmaxdim!(sweeps, 50, 100, 200, 400, 600)
    setcutoff!(sweeps, 1e-12)
    setnoise!(sweeps, 1e-5, 1e-6, 1e-7, 1e-8, 0.0)

    E_0, psi = dmrg(H_mpo, psi_0, sweeps)

    # ============================================
    # Observables: Total Color, Fermion Density, Electric Flux
    # ============================================
    Qx2 = real(inner(Qx, psi, Qx, psi))
    Qy2 = real(inner(Qy, psi, Qy, psi))
    Qz2 = real(inner(Qz, psi, Qz, psi))
    Qtotal2 = Qx2 + Qy2 + Qz2

    println("Ground State Energy E0 = ", E_0)
    println("<Qtotal²> = ", Qtotal2, " (should be ≈ 0)")

    z_values = expect(psi, "Z")
    n_r = [(1.0 - z_values[r(i)]) / 2.0 for i in 0:(N - 1)]
    n_g = [(1.0 - z_values[g(i)]) / 2.0 for i in 0:(N - 1)]
    n_tot = n_r .+ n_g

    flux_squared = zeros(Float64, N - 1)
    for x in 0:(N - 2)
        Lx, Ly, Lz = cumulative_colour_operators(sites, x)
        flux_squared[x + 1] = real(inner(Lx, psi, Lx, psi)) +
                              real(inner(Ly, psi, Ly, psi)) +
                              real(inner(Lz, psi, Lz, psi))
    end

    # ============================================
    # Save Results to HDF5
    # ============================================
    output_file = "SU2_N$(N)_n1$(n1)_n2$(n2).h5"
    h5open(output_file, "w") do f
        write(f, "psi", psi)
        write(f, "E0", E_0)
        write(f, "n_r", n_r)
        write(f, "n_g", n_g)
        write(f, "n_tot", n_tot)
        write(f, "flux_squared", flux_squared)
        write(f, "Qtotal2", Qtotal2)
        write(f, "N", N)
        write(f, "n1", n1)
        write(f, "n2", n2)
    end
    println("Saved results to ", output_file)

    # ============================================
    # Plot Fermion Density & Chromoelectric Field Configuration
    # ============================================
    site_axis = 0:(N - 1)
    link_axis = 0.5:1.0:(N - 1.5)

    p1 = plot(site_axis, n_r;
        marker=:square, label="Red ⟨n_r⟩", color=:crimson,
        ylabel="Occupation", title="Staggered Fermion Occupations (N=$N, m=$m, g=$g_coupling)")
    plot!(p1, site_axis, n_g;
        marker=:circle, label="Green ⟨n_g⟩", color=:seagreen)
    vline!(p1, [n1, n2]; linestyle=:dash, color=:black, label="External charges")

    p2 = plot(site_axis, n_tot;
        marker=:diamond, label="Total ⟨n_r + n_g⟩", color=:navy,
        ylabel="Total Density", title="Local Matter Density per Staggered Site")
    vline!(p2, [n1, n2]; linestyle=:dash, color=:black, label="")

    p3 = plot(link_axis, flux_squared;
        marker=:circle, linewidth=2, label="⟨L²⟩", color=:darkorange,
        xlabel="Lattice Position (staggered index n)", ylabel="⟨L²⟩",
        title="Chromoelectric Flux Profile (String Between n1=$n1 and n2=$n2)")
    vline!(p3, [n1, n2]; linestyle=:dash, color=:black, label="")

    final_plot = plot(p1, p2, p3; layout=(3, 1), size=(850, 950), margin=5Plots.mm)
    savefig(final_plot, "SU2_flux_and_density_N$(N)_n1$(n1)_n2$(n2).png")
    display(final_plot)
end

run_simulation()