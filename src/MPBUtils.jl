module MPBUtils
# ---------------------------------------------------------------------------------------- #

using Crystalline
using Crystalline: TEST_αβγ, DEFAULT_ATOL
using Crystalline: AbstractFourierLattice, 
                   matching_littlegroups, constant, formatirreplabel
using StaticArrays
using Statistics: quantile # for `filling2isoval`
using DocStringExtensions
using Requires
using DelimitedFiles

# ---------------------------------------------------------------------------------------- #

export prepare_mpbcalc,
       prepare_mpbcalc!,
       mpb_calcname,
       write_lgs_to_mpb!,
       lattice_from_mpbparams,
       kvecs_from_mpbparams,
       filling2isoval,
       isoval2filling

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("filling2isoval.jl")
include("export2mpb.jl")
include("read_utils.jl")
export parse_dim, parse_sgnum, load_symdata, symdata2representation, singlebandirreps

# ---------------------------------------------------------------------------------------- #

function __init__()
    # plotting utitilities when PyPlot is loaded
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin  
        include("compat/pyplot.jl")
        export plot_lattice_from_mpbparams
    end
end

# ---------------------------------------------------------------------------------------- #
end