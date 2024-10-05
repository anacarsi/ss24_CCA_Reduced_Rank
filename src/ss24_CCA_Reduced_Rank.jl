# Ana Carsi 2024
module ss24_CCA_Reduced_Rank

using Plots
using Manifolds
using LinearAlgebra
using Manopt
using Random

include("functions.jl")

export cca_on_stiefel, cca_with_tracking, canonical_correlation_heatmap, optimization_path_plot

end
