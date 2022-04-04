module MPBUtils
# ---------------------------------------------------------------------------------------- #

using Crystalline
using Crystalline: AbstractFourierLattice,
                   matching_littlegroups,
                   DEFAULT_ATOL, label, formatirreplabel,
                   TEST_αβγ
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
<<<<<<< HEAD
       isoval2filling,
       extract_candidate_symmetryvectors
=======
       isoval2filling, make_symmetryvectors
>>>>>>> 0be796136771b4b6eb260e854b76754ee28d0029

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("filling2isoval.jl")
include("export2mpb.jl")
include("read_utils.jl")
export read_symdata, extract_multiplicities, extract_individual_multiplicities,
       collect_separable, pick_lgirreps

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