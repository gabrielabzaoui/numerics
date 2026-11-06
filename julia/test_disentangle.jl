using ITensors
using ITensorMPS
using LinearAlgebra
using Serialization

clifford_group=deserialize("two_qubit_clifford_group.jls")
N=18
sites=siteinds("S=1/2",N)
psi=random_mps(sites)
#s1=siteind(psi,1)
#s2=siteind(psi,2)
J=4.0
h=OpSum()

for j in 1:N-1
    global h += J, "Sz", j, "Sz", j+1
    global h += J/2, "S+", j, "S-", j+1
    global h += J/2, "S-", j, "S+", j+1
end

H = MPO(h, sites)
@time psi=psi = tdvp(
        H,
        -1im *0.5,
        psi;
        nsteps = 10,
        maxdim = 200,
        cutoff = 1e-13,
        nsite = 2,
        normalize=true
    )
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
Tensors=[]
psi_opt=copy(psi)
@time for i=3:N 
    psi_iterate=copy(psi)
    s1=siteind(psi_iterate,i-1)#siteind(psi,1)
    s2=siteind(psi_iterate,i)#siteind(psi,2)
    T=nothing
    S_before=entanglement(psi_iterate,i-1)
    S_opt=S_before
    best_state=nothing
    @time for U in clifford_group#[1:720]
        
        U_tensor=ITensor(reshape(U,2,2,2,2),prime(s1),prime(s2),s1,s2)
        #global psi=
        psi_test=apply(U_tensor,copy(psi_opt))
        S_test=entanglement(psi_test,i-1)
        if S_test<S_opt
            best_state = psi_test
            T=U_tensor
            S_opt=S_test
            
            
        end
    end
    if T!==nothing
        push!(Tensors,T)
        global psi_opt=best_state
    end
    println("Bond $i: $S_before -> $S_opt")
end
#psi_opt=copy(psi)
#@time for T in Tensors
 #  global psi_opt=apply(T,psi_opt)
#end
for i=2:N
    display(entanglement(psi,i)-entanglement(psi_opt,i))
    println("\n")
end