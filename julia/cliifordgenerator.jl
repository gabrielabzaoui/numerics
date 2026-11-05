
using LinearAlgebra

# ============================================================
# 1. Define the single-qubit Clifford generators
# ============================================================

I2 = ComplexF64[
    1 0
    0 1
]

H = ComplexF64[
    1  1
    1 -1
] / sqrt(2)

S = ComplexF64[
    1  0
    0 im
]

# ============================================================
# 2. Define the two-qubit generators
# ============================================================

H1 = kron(H, I2)
H2 = kron(I2, H)

S1 = kron(S, I2)
S2 = kron(I2, S)

CNOT = ComplexF64[
    1 0 0 0
    0 1 0 0
    0 0 0 1
    0 0 1 0
]

generators = [H1, H2, S1, S2, CNOT]

# ============================================================
# 3. Remove global phase
# ============================================================

function canonicalize(U; tol=1e-10)
    V = copy(U)

    k = findfirst(x -> abs(x) > tol, V)

    if k === nothing
        error("Zero matrix encountered")
    end

    V .*= conj(V[k]) / abs(V[k])

    # Round to avoid floating-point differences
    V = round.(V; digits=10)

    return V
end

# ============================================================
# 4. Generate the complete two-qubit Clifford group
# ============================================================

function generate_two_qubit_clifford_group(generators)

    Id = Matrix{ComplexF64}(I, 4, 4)

    group = Matrix{ComplexF64}[]
    queue = Matrix{ComplexF64}[]

    seen = Set{Any}()

    Id_canonical = canonicalize(Id)

    push!(group, Id_canonical)
    push!(queue, Id_canonical)
    push!(seen, Tuple(vec(Id_canonical)))

    head = 1

    while head <= length(queue)

        U = queue[head]
        head += 1

        for G in generators

            V = canonicalize(G * U)
            key = Tuple(vec(V))

            if !(key in seen)

                push!(seen, key)
                push!(group, V)
                push!(queue, V)

            end
        end
    end

    return group
end

# ============================================================
# 5. Generate and benchmark
# ============================================================

println("Generating two-qubit Clifford group...")

@time clifford_group =
    generate_two_qubit_clifford_group(generators)

println()
println("Generation completed!")
println("Number of elements: ", length(clifford_group))

# ============================================================
# 6. Validate the group
# ============================================================

println()
println("Checking unitarity...")

max_error = 0.0

for U in clifford_group

    err = norm(U' * U - I)

    global max_error = max(max_error, err)

end

println("Maximum unitarity error: ", max_error)

# ============================================================
# 7. Save the generated group
# ============================================================

using Serialization

serialize("two_qubit_clifford_group.jls", clifford_group)

println()
println("Saved as two_qubit_clifford_group.jls")
