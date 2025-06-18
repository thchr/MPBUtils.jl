# MPBUtils

MPBUtils.jl interfaces with [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to set up and post-process [mpb (MIT Photonic Bands)](https://github.com/NanoComp/mpb) calculations of band connectivity and topology of photonic crystals using symmetry indicators (also known as topological quantum chemistry).

## Installation

This package is not presently registered (and may well change its name in the future). To install it, go to Julia's `pkg>` prompt (by pressing `]`) and type:
```jl
pkg> dev https://github.com/thchr/MPBUtils.jl
```
[SymmetryBases.jl](https://github.com/thchr/SymmetryBases.jl) is a dependency of MPBUtils.jl (which also not registered); by `dev`ing (rather than `add`ing), it will also be automatically installed.

## Functionality

The package at present contains two sets of distinct utilities:

1. Utilities to facilitate the symmetry-based topological analysis of photonic crystals, relying on tools in [Crystalline.jl](https://github.com/thchr/Crystalline.jl) in combination with [mpb](https://github.com/NanoComp/mpb)'s ability to compute the symmetry eigenvalues of the photonic band structure.
2. Exportation and importation of Guile parseable job scripts for [mpb](https://github.com/NanoComp/mpb)'s .ctl interface. This utility is subject to future removal as its effective use requires `.ctl` files that are not included in this repository.

We describe the utilities in point 1 by example below.

## Examples

### 2D photonic crystal

MPBUtils.jl provides a set of convenience tools to initialize and process symmetry analyses of photonic crystal band structures, aimed at making this possible in an interactive manner via [mpb](https://github.com/NanoComp/mpb)'s python interface (called from Julia via PythonCall.jl). To illustrate the functionality, we will first consider a simple 2D photonic crystal example.

First, we make the [mpb](https://github.com/NanoComp/mpb) python interface accessible via Julia and also load the Crystalline.jl and MPBUtils.jl packages:
```jl
# --- load relevant packages ---
using Crystalline, MPBUtils
using PythonCall # use CondaPkg.jl to add meep & mpb (provided by "pymeep") if not already installed
mp = pyimport("meep")
mpb = pyimport("meep.mpb")
```
Note that, in order to compute symmetry eigenvalues via [mpb](https://github.com/NanoComp/mpb)'s  python interface, a relatively recent version of [meep](https://github.com/NanoComp/meep) (≥v1.23.0) is required.

Then we initialize a 2D photonic crystal calculation:
```jl
# --- mpb: geometry & solver initialization ---
ms = mpb.ModeSolver(
    num_bands        = 10,
    resolution       = 32,
    geometry_lattice = mp.Lattice(size=[1,1]),
    geometry         = pylist([mp.Block(center=[0,0], size=[0.3,0.3],  # a 15° rotated square
                                 material=mp.Medium(epsilon=16),
                                 e1=[cosd(15),sind(15)],
                                 e2=[cosd(105),sind(105)])])
    )
ms.init_params(p = mp.TM, reset_fields=true) # solve for TM modes
```

This structure has the symmetry of [plane group 10 (p4)](https://www.cryst.ehu.es/cgi-bin/plane/programs/nph-plane_getgen?gnum=10&type=plane&what=gp). In preparation for the following steps, we first obtain relevant group theory related data for this plane group via Crystalline.jl:
```jl
# --- band representations, littlegroups, & irreps ---
D, sgnum = 2, 10 # dimension and plane group (p4, with Z₂ indicator group)
brs = primitivize(calc_bandreps(sgnum, Val(D))) # elementary band representations
lgirsv = irreps(brs)                            # small irreps & little groups assoc. w/ `brs`
```

We take care above to convert the internally referenced little groups of `brs` to a primitive setting (via `primitivize`), since `calc_bandreps` (and indeed, all other accessors in Crystalline) by default returns operations in a _conventional_ setting; when we consider the symmetry eigenvalues, however, we must work in a _primitive_ setting since our computational unit cell also will be a primitive one. That is, the operations and **k**-points used when calculating symmetry eigenvalues must refer to the same setting as used in the associated unit cell; to avoid redundant band folding, this should be a _primitive_ unit cell.
The reduction step is actually redundant in plane group 10 (p4), as its conventional setting is already primitive. However, we stress it here to emphasize that this is necessary in the general case (for centered Bravais lattices).

Next, using [mpb](https://github.com/NanoComp/mpb), we compute the relevant symmetry eigenvalues of the photonic band structure at each of the **k**-points featured in `brs`, `lgs`, and `lgirsv`:
```jl
# --- compute band symmetry data ---
# symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩, indexed over k-labels `klab`, band indices `n`, and
# operations `gᵢ` (index `i`)
symeigsv = Vector{Vector{Vector{ComplexF64}}}(undef, length(lgirsv)) # symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩
for (kidx, lgirs) in enumerate(lgirsv)
    lg = group(lgirs)
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsv[kidx] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:pyconvert(Int, ms.num_bands)]
    for (i, gᵢ) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(gᵢ))..., [0,0,1]) # decompose gᵢ = {W|w}
        w = mp.Vector3(translation(gᵢ)...)
        symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
        symeigs = pyconvert(Vector{ComplexF64}, symeigs) # convert from Python types to Julia types
        setindex!.(symeigsv[kidx], symeigs, i) # update container of symmetry eigenvalues
    end
end
```

Because the photonic band structure is singular at zero frequency, [mpb](https://github.com/NanoComp/mpb) will not generally be able to assign the appropriate symmetry eigenvalue at (**k** = Γ, ω = 0).
To correct for this, we use MPBUtils.jl's `fixup_gamma_symmetry!` on our symmetry eigenvalue data `symeigsv`:
```jl
# --- fix singular photonic symmetry content at Γ, ω=0 ---
fixup_gamma_symmetry!(symeigsv, lgirsv, :TM) # must specify polarization (:TE or :TM) for `D=2`
```

Finally, we use the elementary band representations and little group irreps to analyze the symmetry eigenvalue data `symeigsv`, extracting the associated band connectivity and band topology of the separable bands in our calculation:
```jl
# --- analyze connectivity and topology of symmetry data ---
summaries = collect_compatible_detailed(symeigsv, brs)
```

For the above structure, this returns the following vector of `BandSummary`s:
```jl
julia> summaries
5-element Vector{BandSummary{2}}:
 [X₁, M₁, Γ₁] (1 band) trivial
 [3X₂, M₂+M₃M₄, Γ₁+Γ₃Γ₄] (3 bands) trivial
 [2X₁, M₃M₄, 2Γ₂] (2 bands) fragile
 [X₁, M₂, Γ₁] (1 band) nontrivial
 [X₁+X₂, M₁+M₂, Γ₃Γ₄] (2 bands) trivial
```

Each band summary contains detailed information about the associated photonic bands. We can e.g., inspect the 4th band grouping in more detail:
```jl
julia> summaries[4]
1-band BandSummary{2}:
 bands:      7:7
 n:          [X₁, M₂, Γ₁] (1 band)
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
As an example, the following scripts sets up a photonic crystal calculation with the symmetry of [space group 81 (P-4)](https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-getgen?what=gp&gnum=81&what=gp) and analyses its band connectivity and topology from symmetry (execution should take on the order of 10-20 seconds):

```jl
# --- load relevant packages ---
using Crystalline, MPBUtils
using PythonCall
mp = pyimport("meep")
mpb = pyimport("meep.mpb")

# --- mpb: geometry & solver initialization ---
r = 0.15
mat = mp.Medium(epsilon=12)
ms = mpb.ModeSolver(
    num_bands        = 20,
    geometry_lattice = mp.Lattice(basis1=[1,0,0], basis2=[0,1,0], basis3=[0,0,1],
                                  basis_size=[1,1,1]),
    geometry         = pylist([
                        mp.Sphere(center=[0.35,0.1,-.2],   radius=r, material=mat),
                        mp.Sphere(center=[-0.1,0.35,.2],   radius=r, material=mat),
                        mp.Sphere(center=[-0.35,-0.1,-.2], radius=r, material=mat),
                        mp.Sphere(center=[0.1,-0.35,.2],   radius=r, material=mat)]),
    resolution       = 16
    )
ms.init_params(p=mp.ALL, reset_fields=true)

# --- band representations, littlegroups, & irreps ---
D, sgnum = 3, 81                                # P-4 (Z₂×Z₂ symmetry indicator group)
brs = primitivize(calc_bandreps(sgnum, Val(D))) # elementary band representations
lgirsv = irreps(brs)                            # associated little groups & small irreps

# --- compute band symmetry data ---
# symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩, indexed over k-labels `klab`, band indices `n`, and
# operations `gᵢ` (index `i`)
symeigsv = Vector{Vector{Vector{ComplexF64}}}(undef, length(lgirsv)) # symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩
for (kidx, lgirs) in enumerate(lgirsv)
    lg = group(lgirs)
    kv = mp.Vector3(position(lg)()...)
    ms.solve_kpoint(kv)

    symeigsv[kidx] = [Vector{ComplexF64}(undef, length(lg)) for n in 1:pyconvert(Int, ms.num_bands)]
    for (i, gᵢ) in enumerate(lg)
        W = mp.Matrix(eachcol(rotation(gᵢ))...) # decompose gᵢ = {W|w}
        w = mp.Vector3(translation(gᵢ)...)
        symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
        symeigs = pyconvert(Vector{ComplexF64}, symeigs) # convert from Python types to Julia types
        setindex!.(symeigsv[kidx], symeigs, i) # update container of symmetry eigenvalues
    end
end

# --- fix singular photonic symmetry content at Γ, ω=0 ---
fixup_gamma_symmetry!(symeigsv, lgirsv)

# --- analyze connectivity and topology of symmetry data ---
summaries = collect_compatible_detailed(symeigsv, brs)
```

Producing the result:
```jl
julia> summaries
5-element Vector{BandSummary{3}}:
 [Z₃Z₄, M₁+M₂, A₁+A₂, X₁+X₂, -Γ₁+Γ₂+Γ₃Γ₄, R₁+R₂] (2 bands) nontrivial
 [Z₃Z₄, M₃M₄, A₃A₄, X₁+X₂, Γ₃Γ₄, R₁+R₂] (2 bands) nontrivial
 [Z₁+Z₂+Z₃Z₄, M₁+M₂+M₃M₄, A₁+A₂+A₃A₄, 2X₁+2X₂, Γ₁+Γ₂+Γ₃Γ₄, 2R₁+2R₂] (4 bands) trivial
 [Z₂+Z₃Z₄, M₁+2M₂, 2A₁+A₂, X₁+2X₂, Γ₂+Γ₃Γ₄, R₁+2R₂] (3 bands) nontrivial
 [Z₁+Z₃Z₄, M₁+M₃M₄, A₂+A₃A₄, 2X₁+X₂, Γ₁+Γ₃Γ₄, 2R₁+R₂] (3 bands) nontrivial
```

## Collaboration and how to cite

Please consider reaching out to us directly if you find the included functionality interesting.
See also the papers below, for which the tools were developed:

- T. Christensen, H.C. Po, J.D. Joannopoulos, & M. Soljačić, *Location and topology of the fundamental gap in photonic crystals*, [Physical Review X **12**, 021066 (2022)](https://doi.org/10.1103/PhysRevX.12.021066).

- A. Ghorashi, S. Vaidya, M.C. Rechtsman, W.A. Benalcazar, M. Soljačić, T. Christensen, *Prevalence of Two-Dimensional Photonic Topology*, [Physical Review Letters **133**, 056602 (2024)](https://doi.org/10.1103/PhysRevLett.133.056602).