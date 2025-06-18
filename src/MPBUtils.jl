module MPBUtils
# ---------------------------------------------------------------------------------------- #

using Crystalline
using Crystalline: AbstractFourierLattice, find_multiplicities # internals
using Crystalline: TEST_αβγ # default value for setting free k-vector parameters `αβγ`

# ---------------------------------------------------------------------------------------- #

using StaticArrays
using DocStringExtensions
using DelimitedFiles
using LinearAlgebra
import SymmetryBases # only used for `calc_detailed_topology`

# ---------------------------------------------------------------------------------------- #

export prepare_mpbcalc
export prepare_mpbcalc!
export mpb_calcname
export write_lgs_to_mpb!
export lattice_from_mpbparams
export kvecs_from_mpbparams

export read_symdata
export extract_multiplicities
export extract_all_multiplicities

export fixup_gamma_symmetry!

export BandSummary
export collect_compatible_detailed

# ---------------------------------------------------------------------------------------- #

@deprecate detailed_symeigs_analysis collect_compatible_detailed

# ---------------------------------------------------------------------------------------- #

include("utils.jl")
include("export2mpb.jl")
include("bandsummary.jl")
include("photonic_symeig_utils.jl")
include("read_utils.jl")

# ---------------------------------------------------------------------------------------- #

end # module MPBUtils