# Ana Carsi 
# CCA with Stiefel manifold optimization
function cca_stiefel(X, Y, k)
    n, p = size(X)
    q = size(Y, 2)
    
    M = Stiefel(p, k) × Stiefel(q, k)
    
    function cost(point)
        U, V = point
        -tr((U' * X' * Y * V) * (U' * X' * Y * V)')
    end
    
    function egrad(point)
        U, V = point
        XY = X' * Y
        gradU = -2 * XY * V * (U' * XY * V)'
        gradV = -2 * XY' * U * (V' * XY' * U)'
        return gradU, gradV
    end
    
    problem = Problem(M, cost, egrad=egrad)
    
    initial_point = (qr(randn(p, k)).Q, qr(randn(q, k)).Q)
    
    result = trust_regions(problem, initial_point)
    
    return result
end
