using ITensors,ITensorMPS
using Plots
using Random
Random.seed!(1234)
begin
    N=4
    sites=siteinds("S=1/2",N)
    h=2
    H=OpSum()
    for i=1:N
        global H+= h,"Sx",i
    
    end
    for i=1:N-1
        global H+=1,"Sz",i,"Sz",i+1
    end
    H+=1,"Sz",N,"Sz",1
    function pauli_tensor(s, P)
    return ITensor(P, prime(s), dag(s))
    end

    const I2 = ComplexF64[1 0; 0 1]
    const X  = ComplexF64[0 1; 1 0]
    const Y  = ComplexF64[0 -im; im 0]
    const Z  = ComplexF64[1 0; 0 -1]

    const PAULIS = [I2, X, Y, Z]
    psi=random_mps(sites)
    normalize!(psi)
    function M_2(psi)
        
        Q=0.0
        for p in Iterators.product(ntuple(i->1:4,N)...)
            P = PAULIS[p[1]]

                for j in 2:N
                    P = kron(P, PAULIS[p[j]])
                end
            Q+= abs2(dot(psi,P*psi))^2

        end
    
        return -log2(Q/2^N)

    end

    function mps_to_vector(psi)

        T = psi[1]

        for j in 2:length(psi)
            T *= psi[j]
        end

        return vec(Array(T, siteinds(psi)...))
    end
    function M_2_mps(psi)
        N = length(psi)
        sites = siteinds(psi)
        E_I=ITensor[]
        E_X=ITensor[]
        E_Y=ITensor[]
        E_Z=ITensor[]
        for i=1:N
            A=psi[i]
            s=sites[i]
            I2_tensor=ITensor(I2,s,prime(s))
            X_tensor=ITensor(X,s,prime(s))
            Y_tensor=ITensor(Y,s,prime(s))
            Z_tensor=ITensor(Z,s,prime(s))
            push!(E_I,dag(A)*I2_tensor*prime(A))
            push!(E_X,dag(A)*X_tensor*prime(A))
            push!(E_Y,dag(A)*Y_tensor*prime(A))
            push!(E_Z,dag(A)*Z_tensor*prime(A))
        end
        
        transfer_tensors=[E_I,E_X,E_Y,E_Z]
        T=[ITensor() for i=1:N]
        
        for i=1:N
            #for p in tuple
                
                T[i]+=E_I[i]*conj(prime(E_I[i],1))*prime(E_I[i],2)*conj(prime(E_I[i],3))
                T[i]+=E_X[i]*conj(prime(E_X[i],1))*prime(E_X[i],2)*conj(prime(E_X[i],3))
                T[i]+=E_Y[i]*conj(prime(E_Y[i],1))*prime(E_Y[i],2)*conj(prime(E_Y[i],3))
                T[i]+=E_Z[i]*conj(prime(E_Z[i],1))*prime(E_Z[i],2)*conj(prime(E_Z[i],3))
            #end
        end
        Q=T[1]
        for i=2:N
            Q*=T[i]
        end
        #println(scalar(Q))
        return -log2(real(scalar(Q))/2^N)
    end
    @time println(M_2(mps_to_vector(psi)))
    @time println(M_2_mps(psi))
end
