module cliffordaug
    
using ITensors
using ITensorMPS
using LinearAlgebra
using Serialization
#sites=siteinds("S=1/2",20)
#psi=MPS(sites,["Up" for i=1:20])
export optimization, entanglement

#X=[0 1.;
 #   1. 0]
#display(X)



I2 = ComplexF64[1 0; 0 1]

H = ComplexF64[1 1; 1 -1] / sqrt(2)

S = ComplexF64[1 0; 0 im]

CNOT = ComplexF64[
    1 0 0 0
    0 1 0 0
    0 0 0 1
    0 0 1 0
]

# Two-qubit generators
H1 = kron(H, I2)
H2 = kron(I2, H)

S1 = kron(S, I2)
S2 = kron(I2, S)
# Controlled-Z
CZ = ComplexF64[
    1 0 0  0
    0 1 0  0
    0 0 1  0
    0 0 0 -1
]

# Controlled-Y (first qubit is control)
CY = ComplexF64[
    1 0 0  0
    0 1 0  0
    0 0 0 -im
    0 0 im  0
]

# Reverse CNOT (second qubit is control)
CNOT_rev = ComplexF64[
    1 0 0 0
    0 0 0 1
    0 0 1 0
    0 1 0 0
]

# SWAP
SWAP = ComplexF64[
    1 0 0 0
    0 0 1 0
    0 1 0 0
    0 0 0 1
]

# iSWAP
iSWAP = ComplexF64[
    1 0 0  0
    0 0 im 0
    0 im 0 0
    0 0 0  1
]

#List = [H1, H2, S1, S2, CNOT,CZ, CY,CNOT_rev,SWAP,iSWAP]
List=deserialize("two_qubit_clifford_group.jls")


function optimization(psi::MPS,N)
    
    for i=1:N
        for i=1:length(psi)-1
            orthogonalize!(psi, i)
            
            theta=psi[i]*psi[i+1]
            s1=siteind(psi,i)
            s2=siteind(psi,i+1)
            left_ind=uniqueinds(psi[i],psi[i+1])
            theta=cliffordaug_sites(psi[i]*psi[i+1],left_ind,s1,s2)
            
            U, S, V=svd(theta,left_ind)    
            psi[i]=U*S
            psi[i+1]=V
        end
        for i in (length(psi)-1):-1:1
            orthogonalize!(psi, i)
            theta=psi[i]*psi[i+1]
            s1=siteind(psi,i)
            s2=siteind(psi,i+1)
            left_ind=uniqueinds(psi[i],psi[i+1])
            theta=cliffordaug_sites(psi[i]*psi[i+1],left_ind,s1,s2)
            
            U, S, V=svd(theta,left_ind)    
            psi[i]=U*S
            psi[i+1]=V
        end
    end    
    return psi
end

function cliffordaug_sites(theta,left_ind,s1,s2)
    S=entropy(theta,left_ind)
    theta_opt=theta
    
    
    
    for U in List
        U_tensor=ITensor(reshape(U,2,2,2,2),prime(s1),prime(s2),s1,s2)
        A=U_tensor*theta
        A = replaceind(A, prime(s1), s1)
        A = replaceind(A, prime(s2), s2)
        S_canidate=entropy(A,left_ind)
        if S_canidate< S
            theta_opt=A
            S=S_canidate
        end
    end
    return theta_opt
end

function entropy(theta,left_inds)
    S_ent=0.
    #left_inds=uniqueinds(psi1,psi2)
    U, S, V=svd(theta,left_inds)
    for n=1:dim(S, 1)
            p = S[n,n]^2
            S_ent -= p * log(p)
        end
    return S_ent    
end

#println("hi")
function entanglement(psi,b)
    orthogonalize!(psi, b)
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end

end