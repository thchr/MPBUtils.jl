module MPBUtilsPyPlotExt
# Plot from MPB-data (::AbstractFourierLattice)
using PyPlot: plot
using MPBUtils: lattice_from_mpbparams

function plot_lattice_from_mpbparams(filepath::String; kws...)
    Rs, flat, isoval, _ = lattice_from_mpbparams(filepath)
    plot(flat, Rs; isoval=isoval, kws...) # via Crystalline.jl extension of PyPlot's `plot`
end

end # module MPBUtilsPyPlotExt