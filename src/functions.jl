# Ana Carsi 2024
# Project: Review on Adaptive Canonical Correlation Analysis
# Description: This file contains the functions to perform Canonical Correlation Analysis (CCA) on the Stiefel Manifold.

using Manifolds, Manopt, LinearAlgebra

#= Define functions to perform Canonical Correlation Analysis (CCA) on the Stiefel Manifold.

Parameters:
----------------
- X: Data matrix of size n x p
- Y: Data matrix of size n x q
- k: Number of canonical components to extract

Functions:
----------------
- cca_on_stiefel: Perform CCA on the Stiefel Manifold
- cca_with_tracking: Perform CCA on the Stiefel Manifold with tracking of objective values
- canonical_correlation_heatmap: Plot the canonical correlation heatmap
- optimization_path_plot: Plot the optimization path on the Stiefel Manifold =#
function cca_on_stiefel(X::Matrix{Float64}, Y::Matrix{Float64}, k::Int)
    n, p = size(X)
    q = size(Y, 2)
    
    M1 = Stiefel(p, k)
    M2 = Stiefel(q, k)
    M = ProductManifold(M1, M2)
    
    # Define the cost function to be minimized
    function cost(M::ProductManifold, point::Tuple{Matrix{Float64}, Matrix{Float64}})
        U, V = point
        cost = -tr((U' * X' * Y * V) * (U' * X' * Y * V)')
        println("The cost function value is $cost")
        return cost
    end

    # Define the Euclidean gradient of the cost function
    function egrad(M, point)
        U, V = point
        XY = X' * Y
        gradU = -2 * XY * V * (U' * XY * V)'
        gradV = -2 * XY' * U * (V' * XY' * U)'
        return (gradU, gradV)
    end

    # Define custom retraction method for the product manifold
    function retract(M, p, X)
        U, V = p
        XU, XV = X
        return (retract(M.manifolds[1], U, XU), retract(M.manifolds[2], V, XV))
    end

    # Define custom vector transport method for the product manifold
    function vector_transport(M, p, X, q)
        U, V = p
        XU, XV = X
        qU, qV = q
        return (vector_transport(M.manifolds[1], U, XU, qU), vector_transport(M.manifolds[2], V, XV, qV))
    end


    initial_point = (rand(M1), rand(M2))
    
    result = trust_regions(
        M,
        cost,
        egrad,
        initial_point;
        retraction_method = retract,
        vector_transport_method = vector_transport
    )
    
    return result.x
end

# Function to track CCA convergence over multiple iterations
#=
function cca_with_tracking(X, Y, k, max_iter=100)
    n, px = size(X)
    n, py = size(Y)
    
    Wx, Wy = cca_on_stiefel(X, Y, k)

    obj_values = Float64[]
    for i in 1:max_iter
        obj = -tr((Wx' * X' * Y * Wy) * (Wx' * X' * Y * Wy)')
        push!(obj_values, obj)
        
        Wx, Wy = cca_on_stiefel(X, Y, k)
    end
    
    plot(1:max_iter, obj_values, xlabel="Iteration", ylabel="Objective Value", 
         title="CCA Convergence on Stiefel Manifold", legend=false)
end

function canonical_correlation_heatmap(X, Y, k)
    Wx, Wy = cca_on_stiefel(X, Y, k)
    
    Ux = X * Wx
    Uy = Y * Wy
    
    corr_matrix = cor(Ux, Uy)
    
    heatmap(corr_matrix, xlabel="X Canonical Variates", ylabel="Y Canonical Variates", 
            title="Canonical Correlations", color=:viridis)
end

function optimization_path_plot(X, Y, k, max_iter=100)
    n, px = size(X)
    n, py = size(Y)
    
    Wx, Wy = cca_on_stiefel(X, Y, k)
    
    path_x = [Wx[1,1]]
    path_y = [Wx[2,1]]
    path_z = [Wx[3,1]]
    
    for i in 1:max_iter
        Wx, Wy = cca_on_stiefel(X, Y, k)
        push!(path_x, Wx[1,1])
        push!(path_y, Wx[2,1])
        push!(path_z, Wx[3,1])
    end
    
    plot3d(path_x, path_y, path_z, xlabel="X", ylabel="Y", zlabel="Z", 
           title="Optimization Path on Stiefel Manifold", legend=false)
end
=#