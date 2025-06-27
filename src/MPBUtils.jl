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

# --- PythonCall init -------------------------------------------------------------------- #

using PythonCall: pynew, pycopy!, pyimport, Py, pylist, pyconvert

const mp = pynew()
const mpb = pynew()
function __init__()
    try # import the mp and mpb libraries
        pycopy!(mp, pyimport("meep"))
        pycopy!(mpb, pyimport("meep.mpb"))
        # unset zero subnormals cf. https://github.com/NanoComp/meep/issues/1708#issuecomment-3008003253
        Base.set_zero_subnormals(false)
    catch
        @warn "mpb or meep could not be imported: associated functionality is nonfunctional"
    end
end

export pylist, pyconvert
export mp, mpb

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
export compute_symmetry_eigenvalues

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