
module cliffordaug

using ITensors
using ITensorMPS
using LinearAlgebra
using Serialization

export optimization, entanglement, precompute_gate_tensors

# Load the Clifford group only once, when the module is loaded.
const CliffordGroup = deserialize(
    joinpath(@__DIR__, "two_qubit_clifford_group.jls")
)

# Check the group size.
if length(CliffordGroup) != 11520
    @warn "Expected 11520 projective two-qubit Clifford elements, got $(length(CliffordGroup)). Check group generation and duplicates."
end


# --------------------------------------------------
# 1. Precompute ITensor representations of the gates
# --------------------------------------------------

function precompute_gate_tensors(psi::MPS)

    N = length(psi)

    gate_cache = Vector{Vector{ITensor}}(undef, N-1)

    for b in 1:N-1

        s1 = siteind(psi, b)
        s2 = siteind(psi, b+1)

        gate_cache[b] = [
            ITensor(
                reshape(G, 2, 2, 2, 2),
                prime(s1), prime(s2), s1, s2
            )
            for G in CliffordGroup
        ]

    end

    return gate_cache
end


# --------------------------------------------------
# 2. Entropy directly from the singular values
# --------------------------------------------------

function entropy_from_singular_values(S)

    n = dim(S, 1)

    norm_squared = sum(
        abs2(S[k,k]) for k in 1:n
    )

    if norm_squared == 0
        error("Zero norm encountered in entropy calculation")
    end

    entropy_value = 0.0

    for k in 1:n

        p = abs2(S[k,k]) / norm_squared

        if p > 0
            entropy_value -= p * log(p)
        end

    end

    return entropy_value
end


# --------------------------------------------------
# 3. Clifford augmentation
# --------------------------------------------------

function cliffordaug_sites(
    theta,
    left_inds,
    s1,
    s2,
    gate_tensors;
    atol=1e-12
)

    # Initial state: no Clifford operation.
    U_best, S_best, V_best = svd(theta, left_inds)

    entropy_best = entropy_from_singular_values(S_best)

    # Test each Clifford operator independently.
    for G_tensor in gate_tensors

        A = G_tensor * theta

        A = replaceind(A, prime(s1), s1)
        A = replaceind(A, prime(s2), s2)

        U, S, V = svd(A, left_inds)

        entropy_candidate = entropy_from_singular_values(S)

        if entropy_candidate < entropy_best - atol

            entropy_best = entropy_candidate

            U_best = U
            S_best = S
            V_best = V

        end

    end

    return U_best, S_best, V_best

end


# --------------------------------------------------
# 4. Optimization sweeps
# --------------------------------------------------

function optimization(
    psi::MPS,
    n_sweeps::Integer,
    gate_cache;
    atol=1e-12
)

    N = length(psi)

    for sweep in 1:n_sweeps

        # Forward sweep
        for b in 1:N-1

            orthogonalize!(psi, b)

            theta = psi[b] * psi[b+1]

            s1 = siteind(psi, b)
            s2 = siteind(psi, b+1)

            left_inds = uniqueinds(psi[b], psi[b+1])

            U, S, V = cliffordaug_sites(
                theta,
                left_inds,
                s1,
                s2,
                gate_cache[b];
                atol=atol
            )

            psi[b] = U * S
            psi[b+1] = V

        end

        # Backward sweep
        for b in (N-1):-1:1

            orthogonalize!(psi, b)

            theta = psi[b] * psi[b+1]

            s1 = siteind(psi, b)
            s2 = siteind(psi, b+1)

            left_inds = uniqueinds(psi[b], psi[b+1])

            U, S, V = cliffordaug_sites(
                theta,
                left_inds,
                s1,
                s2,
                gate_cache[b];
                atol=atol
            )

            psi[b] = U * S
            psi[b+1] = V

        end

    end

    return psi
end


# --------------------------------------------------
# 5. Convenience method
# --------------------------------------------------

function optimization(psi::MPS, n_sweeps::Integer; atol=1e-12)

    gate_cache = precompute_gate_tensors(psi)

    return optimization(
        psi,
        n_sweeps,
        gate_cache;
        atol=atol
    )

end


# --------------------------------------------------
# 6. Entanglement entropy
# --------------------------------------------------

function entanglement(psi::MPS, b::Integer)

    N = length(psi)

    if !(2 <= b <= N)
        throw(ArgumentError("b must be between 2 and $N"))
    end

    orthogonalize!(psi, b)

    left_inds = (linkind(psi, b-1), siteind(psi, b))

    U, S, V = svd(psi[b], left_inds)

    return entropy_from_singular_values(S)

end

end
