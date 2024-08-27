# Ana Carsi 
# CCA with Stiefel manifold optimization
using CCAStiefeldemo
using Random

# Generate sample data
Random.seed!(42)
n, p, q, k = 1000, 50, 40, 5
X = randn(n, p)
Y = randn(n, q)

# Run CCA on Stiefel manifold
result = cca_stiefel(X, Y, k)

println("Optimization result:")
println("Final cost: ", result.cost)
println("Iterations: ", result.iteration)
println("Stopped due to: ", result.stop_reason)

U, V = result.x
println("Shape of U: ", size(U))
println("Shape of V: ", size(V))
