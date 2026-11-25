
using ITensors
using ITensorMPS
using LinearAlgebra

include("magic.jl")

function test_GHZ()

    N = 3
    sites = siteinds("S=1/2", N)

    # Start with |000>
    psi = MPS(sites, ["0", "0", "0"])

    # Hadamard gate on site 1
    H = ComplexF64[
        1  1
        1 -1
    ] / sqrt(2)

    h = ITensor(H, sites[1], prime(sites[1]))

    psi = apply(h, psi; ns=[1])

    # CNOT gates
    psi = apply(op("CNOT", sites[1], sites[2]), psi)

    psi = apply(op("CNOT", sites[2], sites[3]), psi)

    normalize!(psi)

    # Calculate M2
    psi_dense = magic.mps_to_vector(psi)

    M2_dense = magic.M_2(psi_dense, N)
    M2_mps = magic.M_2_mps(psi)

    println("\nGHZ stabilizer test")
    println("---------------------------")

    println("Dense M2 = ", M2_dense)
    println("MPS M2   = ", M2_mps)

    println("Difference = ", abs(M2_dense - M2_mps))

    @assert isapprox(M2_dense, 0.0; atol=1e-10)
    @assert isapprox(M2_mps, 0.0; atol=1e-10)

    println("\nGHZ stabilizer test passed!")

end

test_GHZ()
