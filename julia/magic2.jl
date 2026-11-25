
module magic

using ITensors
using ITensorMPS
using LinearAlgebra

export M_2, M_2_mps, mps_to_vector

const I2 = ComplexF64[1 0; 0 1]
const X  = ComplexF64[0 1; 1 0]
const Y  = ComplexF64[0 -im; im 0]
const Z  = ComplexF64[1 0; 0 -1]

const PAULIS = (I2, X, Y, Z)

"""
    M_2(psi, N=round(Int, log2(length(psi))))

Dense reference implementation of the second stabilizer
Rényi measure.
"""
function M_2(psi, N::Integer=round(Int, log2(length(psi))))

    @assert length(psi) == 2^N

    nrm2 = real(dot(psi, psi))
    @assert nrm2 > 0

    # Normalize explicitly.
    psi = psi / sqrt(nrm2)

    Q = 0.0

    for p in Iterators.product(ntuple(_ -> 1:4, N)...)

        P = PAULIS[p[1]]

        for j in 2:N
            P = kron(P, PAULIS[p[j]])
        end

        expectation = dot(psi, P * psi)

        Q += abs2(expectation)^2
    end

    q = real(Q) / 2^N

    return -log2(q)
end


"""
    mps_to_vector(psi)

Convert an MPS to a dense vector.

Only intended for small-system validation.
"""
function mps_to_vector(psi)

    T = psi[1]

    for j in 2:length(psi)
        T *= psi[j]
    end

    return vec(Array(T, siteinds(psi)...))
end


"""
    M_2_mps(psi; check=true)

Calculate the second stabilizer Rényi measure directly
from an MPS.

The Pauli sum is contracted locally using four replicas.

No dense state vector is constructed.

The input MPS may be unnormalized. Its normalization is
accounted for explicitly.
"""
function M_2_mps(psi; check=true)

    N = length(psi)

    @assert N >= 1

    sites = siteinds(psi)

    # --------------------------------------------------
    # 1. Calculate the MPS norm.
    # --------------------------------------------------

    nrm = norm(psi)

    if !(isfinite(nrm) && nrm > 0)
        error("The MPS must have a finite, nonzero norm.")
    end

    # --------------------------------------------------
    # 2. Construct the local four-replica transfer tensors.
    # --------------------------------------------------

    W = Vector{ITensor}(undef, N)

    for i in 1:N

        A = psi[i]
        s = sites[i]

        # E[p] = A† P A
        #
        # Physical indices are contracted locally.
        # Virtual indices remain open.

        E = Vector{ITensor}(undef, 4)

        for p in 1:4

            op = ITensor(PAULIS[p], s, prime(s))

            E[p] = dag(A) * op * prime(A)

        end

        # Four replicas:
        #
        #   E
        #   conj(E), shifted by 2
        #   E,       shifted by 4
        #   conj(E), shifted by 6
        #
        # Every replica has independent virtual indices.

        W[i] = (
            E[1] *
            dag(prime(E[1], 2)) *
            prime(E[1], 4) *
            dag(prime(E[1], 6))
        )

        for p in 2:4

            W[i] += (
                E[p] *
                dag(prime(E[p], 2)) *
                prime(E[p], 4) *
                dag(prime(E[p], 6))
            )

        end

    end

    # --------------------------------------------------
    # 3. Contract the four-replica MPS network.
    # --------------------------------------------------

    T = W[1]

    for i in 2:N

        T *= W[i]

    end

    # An open-boundary MPS must produce a scalar.
    if !isempty(inds(T))
        error(
            "Four-replica contraction left open indices: $(inds(T))"
        )
    end

    Q_raw = real(scalar(T))

    # --------------------------------------------------
    # 4. Correct for possible input normalization.
    #
    # Each expectation value scales as nrm^2.
    # Its fourth power therefore scales as nrm^8.
    # --------------------------------------------------

    Q = Q_raw / nrm^8

    q = Q / 2.0^N

    # --------------------------------------------------
    # 5. Numerical consistency checks.
    # --------------------------------------------------

    if check

        if !isfinite(q) || q <= 0
            error("Invalid stabilizer Rényi argument q = $q")
        end

        tol = 1e-10

        if q > 1 + tol
            error(
                "q = $q exceeds 1. Check the MPS contraction, " *
                "site indices, and normalization."
            )
        end

    end

    # Protect against tiny floating-point overshoots.
    q = min(q, 1.0)

    return -log2(q)

end

end
