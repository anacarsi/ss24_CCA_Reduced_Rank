# Ana Carsi 
# CCA with Stiefel manifold optimization
# scripts/run_cca.jl


using Pkg
Pkg.activate(dirname(@__DIR__)) # Activate the environment of the project's root directory

using ss24_CCA_Reduced_Rank
using Random
using Plots
using Manifolds
using LinearAlgebra
using Manopt

Random.seed!(42)

X = randn(100, 10)
Y = randn(100, 15)
k = 3

p = cca_on_stiefel(X, Y, k)
display(p)
savefig(p, "cca_convergence.png")
