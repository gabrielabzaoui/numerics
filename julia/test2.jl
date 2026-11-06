
using LinearAlgebra
using Serialization

# Load complete Clifford group
clifford_group = deserialize("two_qubit_clifford_group.jls")

# Pauli matrices
I2 = ComplexF64[1 0; 0 1]
X  = ComplexF64[0 1; 1 0]
Y  = ComplexF64[0 -im; im 0]
Z  = ComplexF64[1 0; 0 -1]

# Generate all 16 two-qubit Pauli matrices
paulis = [
    kron(P1, P2)
    for P1 in (I2, X, Y, Z), P2 in (I2, X, Y, Z)
]
paulis = vec(paulis)


function canonical_key(M; tol=1e-8)
    v = vec(ComplexF64.(M))

    idx = findfirst(z -> abs(z) > tol, v)
    isnothing(idx) && error("Zero matrix encountered")

    # Remove global phase
    v ./= v[idx] / abs(v[idx])

    # Convert complex entries into comparable real tuples
    return Tuple(
        Iterators.flatten(
            ((round(real(z), digits=8),
              round(imag(z), digits=8)) for z in v)
        )
    )
end

function coset_key(C, paulis)
    keys = [canonical_key(P * C) for P in paulis]
    return minimum(keys)
end

# Keep one representative per class
clifford_720 = Matrix{ComplexF64}[]
seen = Set{Any}()

for C in clifford_group
    key = coset_key(C, paulis)

    if !(key in seen)
        push!(seen, key)
        push!(clifford_720, ComplexF64.(C))
    end
end

println("Original gates: ", length(clifford_group))
println("Filtered gates: ", length(clifford_720))

# Save
serialize("two_qubit_clifford_720.jls", clifford_720)

@assert length(clifford_720) == 720
