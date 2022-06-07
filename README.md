# MPBUtils

MPBUtils.jl interfaces with [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to set up and post-process [mpb (MIT Photonic Bands)](https://github.com/NanoComp/mpb) calculations of band connectivity and topology of photonic crystals using symmetry indicators (also known as topological quantum chemistry).

## Installation

The package is not presently registered (and may well change its name in the future). To install it, go to Julia's `pkg>` prompt (by pressing `]`) and type:
```jl
pkg> add https://github.com/thchr/MPBUtils.jl
```

## Functionality

The package at present contains two sets of distinct utilities:

1. Utilities to perform band symmetry analysis of photonic structures, assuming the ability to compute the symmetry eigenvalues of the associated photonic band structure (MPB e.g. has this capability).
2. Exportation and importation of Guile parseable job scripts for [mpb's](https://github.com/NanoComp/mpb) .ctl interface. This utility is subject to future removal as its effective use requires `.ctl` files that are not included in this repository.

We describe the utilities in point 1 by example below.

## Examples

### 2D photonic crystal

MPBUtils.jl provides a set of convenience tools to initialize and process symmetry analyses of photonic crystal band structures, aimed at making this possible in an interactive manner via [mpb's](https://github.com/NanoComp/mpb) python interface (called from Julia via PyCall.jl). To illustrate the functionality, we will first consider a simple 2D photonic crystal example.

First, we make the [mpb](https://github.com/NanoComp/mpb) python interface accessible via Julia and also load the Crystalline.jl and MPBUtils.jl packages:
```jl
# --- load relevant packages ---
using Crystalline, MPBUtils
using PyCall
mp = pyimport("meep")
mpb = pyimport("meep.mpb")
```
Note that, in order to compute symmetry eigenvalues via [mpb's](https://github.com/NanoComp/mpb)  python interface, a relatively recent version of [meep](https://github.com/NanoComp/meep) (≥v1.23.0) is required.

Then we initialize a 2D photonic crystal calculation:
```jl
# --- mpb: geometry & solver initialization ---
ms = mpb.ModeSolver(
        num_bands        = 10,
        resolution       = 32,
        geometry_lattice = mp.Lattice(size=[1,1]),
        geometry         = [mp.Block(center=[0,0], size=[0.3,0.3],  # a 15° rotated square
                                     material=mp.Medium(epsilon=16),
                                     e1 = [cosd(15), sind(15)],
                                     e2 = [cosd(105), sind(105)])]
        )
ms.init_params(p = mp.TM, reset_fields=true) # solve for TM modes
```

This structure has the symmetry of plane group p4. In preparation for the following steps, we request relevant group theory related data structures for this plane group via Crystalline.jl:
```jl
# --- band representations, littlegroups, & irreps ---
D, sgnum = 2, 10 # dimension and plane group (p4, with Z₂ indicator group)
brs = bandreps(sgnum, D)                          # elementary band representations
lgs = littlegroups(sgnum, Val(D))                 # little groups
filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
lgirsd = pick_lgirreps(lgs; timereversal=true)    # small irreps associated with `lgs`
```

Next, using [mpb](https://github.com/NanoComp/mpb), we compute the relevant symmetry eigenvalues of the photonic band structure at each of the **k**-points featured in `brs`, `lgs`, and `lgirsd`:
```jl
# --- compute band symmetry data ---
# symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩, indexed over k-labels `klab`, band indices `n`, and
# operations `gᵢ` (index `i`)
symeigsd = Dict{String, Vector{Vector{ComplexF64}}}() # symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩
for (klab, lg) in lgs
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
    for (i, gᵢ) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(gᵢ))..., [0,0,1]) # decompose gᵢ = {W|w}
        w = mp.Vector3(translation(gᵢ)...)
        symeigs = ms.compute_symmetries(W, w)  # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
        setindex!.(symeigsd[klab], symeigs, i) # update container of symmetry eigenvalues
    end
end
```

Because the photonic band structure is singular at zero frequency, [mpb](https://github.com/NanoComp/mpb) will not generally be able to assign the appropriate symmetry eigenvalue at (**k** = Γ, ω = 0).
To correct for this, we use MPBUtils.jl's `fixup_gamma_symmetry!` on our symmetry eigenvalue data `symeigsd`:
```jl
# --- fix singular photonic symmetry content at Γ, ω=0 ---
fixup_gamma_symmetry!(symeigsd, lgs, :TM) # must specify polarization (:TE or :TM) for `D=2`
```

Finally, we use the elementary band representations and little group irreps to analyze the symmetry eigenvalue data `symeigsd`, extracting the associated band connectivity and band topology of the separable bands in our calculation:
```jl
# --- analyze connectivity and topology of symmetry data ---
summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)
```

For the above structure, this returns the following vector of `BandSummaries`:
```jl
julia> summaries
5-element Vector{BandSummary}:
 1-band (trivial): [X₁, M₁, Γ₁]
 3-band (trivial): [3X₂, M₂+M₃M₄, Γ₁+Γ₃Γ₄]
 2-band (fragile): [2X₁, M₃M₄, 2Γ₂]
 1-band (nontrivial): [X₁, M₂, Γ₁]
 2-band (trivial): [X₁+X₂, M₁+M₂, Γ₃Γ₄]
```

Each band summary contains detailed information about the associated photonic bands. We can e.g., inspect the 4th band grouping in more detail:
```jl
julia> summaries[4]
1-band BandSummary:
 bands:      7:7
 n:          X₁, M₂, Γ₁
 topology:   nontrivial
 indicators: 1 ∈ Z₂
```

Adjacent bands can be "stacked" by addition. E.g., to evaluate the topology of the first three band groupings, we can evaluate:
```jl
julia> summaries[1] + summaries[2] + summaries[3] # or simply, sum(summaries[1:3])
6-band BandSummary:
 bands:      1:6
 n:          3X₁+3X₂, M₁+M₂+2M₃M₄, 2Γ₁+2Γ₂+Γ₃Γ₄
 topology:   trivial
```
From which we see that the fragile bands in the 3rd band grouping are trivialized by the trivial bands in the 1st and 2nd band groupings.

### 3D photonic crystal

Analysis of 3D photonic crystals proceeds similarly.
As an example, the following scripts sets up a photonic crystal calculation with the symmetry of space group 81 (P-4) and analyses its band connectivity and topology from symmetry (execution should take on the order of 10-20 seconds):

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
sgnum = 81                                        # P-4 (Z₂×Z₂ symmetry indicator group)
brs = bandreps(sgnum)                             # elementary band representations
lgs = littlegroups(sgnum)                         # little groups
filter!(((klab, _),) -> klab ∈ klabels(brs), lgs) # restrict to k-points in `brs`
lgirsd = pick_lgirreps(lgs; timereversal=true)    # small irreps associated with `lgs`

# --- compute band symmetry data ---
# symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩, indexed over k-labels `klab`, band indices `n`, and
# operations `gᵢ` (index `i`)
symeigsd = Dict{String, Vector{Vector{ComplexF64}}}()
for (klab, lg) in lgs
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsd[klab] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:ms.num_bands]
    for (i, gᵢ) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
        w = mp.Vector3(translation(gᵢ)...)
        symeigs = ms.compute_symmetries(W, w)   # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
        setindex!.(symeigsd[klab], symeigs, i)  # update container of symmetry eigenvalues
    end
end

# --- fix singular photonic symmetry content at Γ, ω=0 ---
fixup_gamma_symmetry!(symeigsd, lgs)

# --- analyze connectivity and topology of symmetry data ---
summaries = analyze_symmetry_data(symeigsd, lgirsd, brs)`
```

Producing the result:
```jl
julia> summaries
5-element Vector{BandSummary}:
 2-band (nontrivial): [A₁+A₂, -Γ₁+Γ₂+Γ₃Γ₄, M₁+M₂, Z₃Z₄, R₁+R₂, X₁+X₂]
 2-band (nontrivial): [A₃A₄, Γ₃Γ₄, M₃M₄, Z₃Z₄, R₁+R₂, X₁+X₂]
 4-band (trivial): [A₁+A₂+A₃A₄, Γ₁+Γ₂+Γ₃Γ₄, M₁+M₂+M₃M₄, Z₁+Z₂+Z₃Z₄, 2R₁+2R₂, 2X₁+2X₂]
 3-band (nontrivial): [2A₁+A₂, Γ₂+Γ₃Γ₄, M₁+2M₂, Z₂+Z₃Z₄, R₁+2R₂, X₁+2X₂]
 3-band (nontrivial): [A₂+A₃A₄, Γ₁+Γ₃Γ₄, M₁+M₃M₄, Z₁+Z₃Z₄, 2R₁+R₂, 2X₁+X₂]
```

## Collaboration and how to cite

Since this package is still in active development, please consider reaching out to us directly if you find the included functionality interesting.
See also the paper below, for which some of the included functionality was developed:

- T. Christensen, H.C. Po, J.D. Joannopoulos, & M. Soljačić, *Location and topology of the fundamental gap in photonic crystals,* [arXiv:2106.10267 (2021)](https://arxiv.org/abs/2106.10267)