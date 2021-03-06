# MPBUtils

MPBUtils.jl interfaces with [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to set up and post-process [mpb (MIT Photonic Bands)](https://github.com/NanoComp/mpb) calculations of band connectivity and topology of photonic crystals using symmetry indicators (also known as topological quantum chemistry).

## Installation

The package is not presently registered (and may well change its name in the future). To install it, go to Julia's `pkg>` prompt (by pressing `]`) and type:
```jl
pkg> add https://github.com/thchr/SymmetryBases.jl
pkg> add https://github.com/thchr/MPBUtils.jl
```
[SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl) is a dependency of MPBUtils.jl (which also not registered, and so also requires manual installation).

## Functionality

The package at present contains two sets of distinct utilities:

1. Utilities to perform band symmetry analysis of photonic structures, assuming the ability to compute the symmetry eigenvalues of the associated photonic band structure (MPB e.g. has this capability).
2. Exportation and importation of Guile parseable job scripts for [mpb](https://github.com/NanoComp/mpb)'s .ctl interface. This utility is subject to future removal as its effective use requires `.ctl` files that are not included in this repository.

We describe the utilities in point 1 by example below.

## Examples

### 2D photonic crystal

MPBUtils.jl provides a set of convenience tools to initialize and process symmetry analyses of photonic crystal band structures, aimed at making this possible in an interactive manner via [mpb](https://github.com/NanoComp/mpb)'s python interface (called from Julia via PyCall.jl). To illustrate the functionality, we will first consider a simple 2D photonic crystal example.

First, we make the [mpb](https://github.com/NanoComp/mpb) python interface accessible via Julia and also load the Crystalline.jl and MPBUtils.jl packages:
```jl
# --- load relevant packages ---
using Crystalline, MPBUtils
using PyCall
mp = pyimport("meep")
mpb = pyimport("meep.mpb")
```
Note that, in order to compute symmetry eigenvalues via [mpb](https://github.com/NanoComp/mpb)'s  python interface, a relatively recent version of [meep](https://github.com/NanoComp/meep) (???v1.23.0) is required.

Then we initialize a 2D photonic crystal calculation:
```jl
# --- mpb: geometry & solver initialization ---
ms = mpb.ModeSolver(
        num_bands        = 10,
        resolution       = 32,
        geometry_lattice = mp.Lattice(size=[1,1]),
        geometry         = [mp.Block(center=[0,0], size=[0.3,0.3],  # a 15?? rotated square
                                     material=mp.Medium(epsilon=16),
                                     e1 = [cosd(15), sind(15)],
                                     e2 = [cosd(105), sind(105)])]
        )
ms.init_params(p = mp.TM, reset_fields=true) # solve for TM modes
```

This structure has the symmetry of [plane group 10 (p4)](https://www.cryst.ehu.es/cgi-bin/plane/programs/nph-plane_getgen?gnum=10&type=plane&what=gp). In preparation for the following steps, we first obtain relevant group theory related data for this plane group via Crystalline.jl:
```jl
# --- band representations, littlegroups, & irreps ---
D, sgnum = 2, 10 # dimension and plane group (p4, with Z??? indicator group)
brs = bandreps(sgnum, D)                          # elementary band representations
lgs = littlegroups(sgnum, Val(D))                 # little groups
filter!(((klab, _),) -> klab ??? klabels(brs), lgs) # restrict to k-points in `brs`
map!(primitivize, values(lgs))                    # convert to primitive setting
lgirsd = pick_lgirreps(lgs; timereversal=true)    # small irreps associated with `lgs`
```

Note that we convert the little group operations and **k**-points in `lgs` from a conventional to a primitive setting via `primitivize`; this is redundant in plane group 10 (p4), as its conventional setting already primitive. We note it explicitly here to emphasize that this is necessary in the general case (for centered Bravais lattices): symmetry eigenvalues should be computed for a primitive unit cell, with operations and **k**-points referred to the corresponding primitive basis.

Next, using [mpb](https://github.com/NanoComp/mpb), we compute the relevant symmetry eigenvalues of the photonic band structure at each of the **k**-points featured in `brs`, `lgs`, and `lgirsd`:
```jl
# --- compute band symmetry data ---
# symmetry eigenvalues ???E??????|g???D?????????, indexed over k-labels `klab`, band indices `n`, and
# operations `g???` (index `i`)
symeigsd = Dict{String, Vector{Vector{ComplexF64}}}() # symmetry eigenvalues ???E??????|g???D?????????
for (klab, lg) in lgs
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
    for (i, g???) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(g???))..., [0,0,1]) # decompose g??? = {W|w}
        w = mp.Vector3(translation(g???)...)
        symeigs = ms.compute_symmetries(W, w)  # compute ???E??????|g???D????????? for all bands
        setindex!.(symeigsd[klab], symeigs, i) # update container of symmetry eigenvalues
    end
end
```

Because the photonic band structure is singular at zero frequency, [mpb](https://github.com/NanoComp/mpb) will not generally be able to assign the appropriate symmetry eigenvalue at (**k** = ??, ?? = 0).
To correct for this, we use MPBUtils.jl's `fixup_gamma_symmetry!` on our symmetry eigenvalue data `symeigsd`:
```jl
# --- fix singular photonic symmetry content at ??, ??=0 ---
fixup_gamma_symmetry!(symeigsd, lgs, :TM) # must specify polarization (:TE or :TM) for `D=2`
```

Finally, we use the elementary band representations and little group irreps to analyze the symmetry eigenvalue data `symeigsd`, extracting the associated band connectivity and band topology of the separable bands in our calculation:
```jl
# --- analyze connectivity and topology of symmetry data ---
summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)
```

For the above structure, this returns the following vector of `BandSummary`s:
```jl
julia> summaries
5-element Vector{BandSummary}:
 1-band (trivial): [X???, M???, ?????]
 3-band (trivial): [3X???, M???+M???M???, ?????+??????????]
 2-band (fragile): [2X???, M???M???, 2?????]
 1-band (nontrivial): [X???, M???, ?????]
 2-band (trivial): [X???+X???, M???+M???, ??????????]
```

Each band summary contains detailed information about the associated photonic bands. We can e.g., inspect the 4th band grouping in more detail:
```jl
julia> summaries[4]
1-band BandSummary:
 bands:      7:7
 n:          X???, M???, ?????
 topology:   nontrivial
 indicators: 1 ??? Z???
```

Adjacent bands can be "stacked" by addition. E.g., to evaluate the topology of the first three band groupings, we can evaluate:
```jl
julia> summaries[1] + summaries[2] + summaries[3] # or simply, sum(summaries[1:3])
6-band BandSummary:
 bands:      1:6
 n:          3X???+3X???, M???+M???+2M???M???, 2?????+2?????+??????????
 topology:   trivial
```
From which we see that the fragile bands in the 3rd band grouping are trivialized by the trivial bands in the 1st and 2nd band groupings.

### 3D photonic crystal

Analysis of 3D photonic crystals proceeds similarly.
As an example, the following scripts sets up a photonic crystal calculation with the symmetry of [space group 81 (P-4)](https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?what=gp&gnum=81&what=gp) and analyses its band connectivity and topology from symmetry (execution should take on the order of 10-20 seconds):

```jl
# --- load relevant packages ---
using Crystalline, MPBUtils
using PyCall
mp = pyimport("meep")
mpb = pyimport("meep.mpb")

# --- mpb: geometry & solver initialization ---
r = 0.15
mat = mp.Medium(epsilon=12)
ms = mpb.ModeSolver(
    num_bands        = 20,
    geometry_lattice = mp.Lattice(basis1=[1,0,0], basis2=[0,1,0], basis3=[0,0,1],
                                  basis_size=[1,1,1]),
    geometry         = [mp.Sphere(center=[0.35,0.1,-.2],   radius=r, material=mat),
                        mp.Sphere(center=[-0.1,0.35,.2],   radius=r, material=mat),
                        mp.Sphere(center=[-0.35,-0.1,-.2], radius=r, material=mat),
                        mp.Sphere(center=[0.1,-0.35,.2],   radius=r, material=mat)],
    resolution       = 16
    )
ms.init_params(p=mp.ALL, reset_fields=true)

# --- band representations, littlegroups, & irreps ---
sgnum = 81                                        # P-4 (Z?????Z??? symmetry indicator group)
brs = bandreps(sgnum)                             # elementary band representations
lgs = littlegroups(sgnum)                         # little groups
filter!(((klab, _),) -> klab ??? klabels(brs), lgs) # restrict to k-points in `brs`
map!(primitivize, values(lgs))                    # convert to primitive setting
lgirsd = pick_lgirreps(lgs; timereversal=true)    # small irreps associated with `lgs`

# --- compute band symmetry data ---
# symmetry eigenvalues ???E??????|g???D?????????, indexed over k-labels `klab`, band indices `n`, and
# operations `g???` (index `i`)
symeigsd = Dict{String, Vector{Vector{ComplexF64}}}()
for (klab, lg) in lgs
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
    for (i, g???) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(g???))...) # decompose g??? = {W|w}
        w = mp.Vector3(translation(g???)...)
        symeigs = ms.compute_symmetries(W, w)   # compute ???E??????|g???D????????? for all bands
        setindex!.(symeigsd[klab], symeigs, i)  # update container of symmetry eigenvalues
    end
end

# --- fix singular photonic symmetry content at ??, ??=0 ---
fixup_gamma_symmetry!(symeigsd, lgs)

# --- analyze connectivity and topology of symmetry data ---
summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)
```

Producing the result:
```jl
julia> summaries
5-element Vector{BandSummary}:
 2-band (nontrivial): [A???+A???, -?????+?????+??????????, M???+M???, Z???Z???, R???+R???, X???+X???]
 2-band (nontrivial): [A???A???, ??????????, M???M???, Z???Z???, R???+R???, X???+X???]
 4-band (trivial): [A???+A???+A???A???, ?????+?????+??????????, M???+M???+M???M???, Z???+Z???+Z???Z???, 2R???+2R???, 2X???+2X???]
 3-band (nontrivial): [2A???+A???, ?????+??????????, M???+2M???, Z???+Z???Z???, R???+2R???, X???+2X???]
 3-band (nontrivial): [A???+A???A???, ?????+??????????, M???+M???M???, Z???+Z???Z???, 2R???+R???, 2X???+X???]
```

## Collaboration and how to cite

Since this package is still in active development, please consider reaching out to us directly if you find the included functionality interesting.
See also the paper below, for which some of the included functionality was developed:

- T. Christensen, H.C. Po, J.D. Joannopoulos, & M. Solja??i??, *Location and topology of the fundamental gap in photonic crystals*, [arXiv:2106.10267 (2021)](https://arxiv.org/abs/2106.10267)