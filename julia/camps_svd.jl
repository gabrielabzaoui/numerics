using ITensors
using ITensorMPS
using LinearAlgebra
using Serialization

clifford_group=deserialize("two_qubit_clifford_720.jls")
N=24
sites=siteinds("S=1/2",N)
state=[isodd(i) ? "Up" : "Dn" for i=1:N]
psi=MPS(sites,state)#random_mps(sites)
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
@time psi=tdvp(
        H,
        -1im *0.5,
        psi;
        nsteps = 100,
        maxdim = 500,
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
function entropy(S)
    vals = [abs2(S[n,n]) for n in 1:dim(S,1)]
    p = vals ./ sum(vals)

    return -sum(x -> x > 0 ? x * log(x) : 0.0, p)
end
function camps(psi::MPS;cutoff=0,sweeps=1)

    clifford_tensors=[]
    for a=1:sweeps    
        for i=1:length(psi)-1
            orthogonalize!(psi, i)
            left_inds=uniqueinds(psi[i],psi[i+1])
            s1=siteind(psi,i)
            s2=siteind(psi,i+1)
            theta=psi[i]*psi[i+1]
            T=nothing
            theta_test=nothing
            U_opt,S_opt,V_opt=svd(theta,left_inds,cutoff=cutoff)
            entropy_opt=entropy(S_opt)
            for U in clifford_group
            U_tensor=ITensor(reshape(U,2,2,2,2),prime(s1),prime(s2),s1,s2)
            
                theta_test=apply(U_tensor,theta) 
                U,S,V=svd(theta_test,left_inds,cutoff=cutoff)
                entropy2=entropy(S)
                if entropy2<entropy_opt
                    U_opt=U
                    S_opt=S
                    V_opt=V
                    T=U_tensor
                    
                end
            end
            if T!==nothing
                push!(clifford_tensors,T)
                psi[i]=U_opt*S_opt
                psi[i+1]=V_opt
            end
        end
        for i=length(psi)-1:-1:1
            orthogonalize!(psi, i)
            left_inds=uniqueinds(psi[i],psi[i+1])
            s1=siteind(psi,i)
            s2=siteind(psi,i+1)
            theta=psi[i]*psi[i+1]
            T=nothing
            theta_test=nothing
            U_opt,S_opt,V_opt=svd(theta,left_inds,cutoff=cutoff)
            entropy_opt=entropy(S_opt)
            for U in clifford_group
            U_tensor=ITensor(reshape(U,2,2,2,2),prime(s1),prime(s2),s1,s2)
            
                theta_test=apply(U_tensor,theta) 
                U,S,V=svd(theta_test,left_inds,cutoff=cutoff)
                entropy2=entropy(S)
                if entropy2<entropy_opt
                    U_opt=U
                    S_opt=S
                    V_opt=V
                    T=U_tensor
                    
                end
            end
            if T!==nothing
                push!(clifford_tensors,T)
                psi[i]=U_opt*S_opt
                psi[i+1]=V_opt
            end
        end
    end
    return clifford_tensors, psi

end
psi_opt=copy(psi)
@time tensors,psi_opt=camps(psi_opt;cutoff=0,sweeps=3)
for i=2:N
display(entanglement(psi,i)-entanglement(psi_opt,i))
end
