# Plot from MPB-data (::AbstractFourierLattice) -------------------------------------------
using .PyPlot

function plot_lattice_from_mpbparams(filepath::String; kwargs...)
    Rs, flat, isoval, _ = lattice_from_mpbparams(filepath)
    plot(flat, Rs; isoval=isoval, kwargs...) # via Crystalline.jl overload of PyPlot.jl's `plot`
end