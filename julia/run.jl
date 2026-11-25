include("su2again.jl")

using ITensors, ITensorMPS
using Plots

N = 6
n1 = 2
n2 = 4

m = 1.9
gauge_g = 2.0
a = 1.0
epsilon = 1.0
lambda = 30.0

energy, psi, n_r, n_g, flux =
    Start_simulation(epsilon,m,gauge_g,a,N,n1,n2,lambda)