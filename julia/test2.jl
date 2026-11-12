
using LinearAlgebra
using Serialization

function generate_camps_cliffords()

    # --------------------------------------------------
    # 1. Single-qubit Clifford generators
    # --------------------------------------------------

    I2 = Matrix{ComplexF64}(I, 2, 2)

    H = ComplexF64[
        1  1
        1 -1
    ] / sqrt(2)

    S = ComplexF64[
        1 0
        0 im
    ]

    # Two-qubit generators
    H1 = kron(H, I2)
    H2 = kron(I2, H)

    S1 = kron(S, I2)
    S2 = kron(I2, S)

    # CNOT: first qubit controls second
    CNOT12 = ComplexF64[
        1 0 0 0
        0 1 0 0
        0 0 0 1
        0 0 1 0
    ]

    # CNOT: second qubit controls first
    SWAP = ComplexF64[
        1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1
    ]

    CNOT21 = SWAP * CNOT12 * SWAP

    generators = [H1, H2, S1, S2, CNOT12, CNOT21]

    # --------------------------------------------------
    # 2. Construct the 16 two-qubit Pauli operators
    # --------------------------------------------------

    X = ComplexF64[
        0 1
        1 0
    ]

    Y = ComplexF64[
        0 -im
        im 0
    ]

    Z = ComplexF64[
        1 0
        0 -1
    ]

    single_paulis = [I2, X, Y, Z]

    paulis = [
        kron(P1, P2)
        for P1 in single_paulis, P2 in single_paulis
    ]

    paulis = vec(paulis)

    # --------------------------------------------------
    # 3. Identify a Clifford modulo Pauli operators
    # --------------------------------------------------

    function clifford_key(C)

        key = Int[]

        for P in paulis

            # Clifford conjugation
            Q = C * P * C'

            # Identify which Pauli operator Q corresponds to.
            # Ignore its sign.
            overlaps = [
                abs(tr(P2' * Q))
                for P2 in paulis
            ]

            push!(key, argmax(overlaps))
        end

        return Tuple(key)
    end

    # --------------------------------------------------
    # 4. Breadth-first generation of the quotient group
    # --------------------------------------------------

    identity_gate = Matrix{ComplexF64}(I, 4, 4)

    representatives = Matrix{ComplexF64}[]
    queue = Matrix{ComplexF64}[]

    seen = Set{Any}()

    key = clifford_key(identity_gate)

    push!(seen, key)
    push!(queue, identity_gate)

    head = 1

    while head <= length(queue)

        C = queue[head]
        head += 1

        for G in generators

            Cnew = G * C

            key = clifford_key(Cnew)

            if !(key in seen)

                push!(seen, key)
                push!(queue, Cnew)

            end
        end

    end

    representatives = queue

    println("Generated Clifford classes: ", length(representatives))

    return representatives
end


# --------------------------------------------------
# 5. Generate and save
# --------------------------------------------------

clifford_720 = generate_camps_cliffords()

@assert length(clifford_720) == 720

serialize("two_qubit_clifford_720.jls", clifford_720)

println("Successfully saved 720 CAMPS Clifford representatives!")
