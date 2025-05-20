using DelimitedFiles

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return the symmetry eigenvalues and little groups as **k**-label indexed `Dict`s for
an MPB symmetry calculation with ID `calcname`.

## Keyword arguments
- `sgnum`, `dir`: If unspecified (or `nothing`), the space group number and dimension must
  be inferrable from `calcname` (`calcname` must then follow the conventions of
  [`mpb_calcname`](@ref)).
- `dir` (default, `"."`): directory location where `calcname`-*.out files can be found.
- `αβγ` (default, `$TEST_αβγ`): free-parameters used for **k**-vectors in setting up the
  MPB calculation in [`prepare_mpbcalc`](@ref).
- `isprimitive` (default, `true`): whether the calculation is in a primitive setting.
- `flip_ksign` (default, `false`): flip the sign of the **k**-vector used in the MPB 
  calculation; can be necessary to match phase conventions in Crystalline.jl.
"""
function read_symdata(
    calcname::AbstractString; 
    sgnum::Union{Int, Nothing} = nothing,
    D::Union{Int, Nothing} = nothing,
    dir::AbstractString = ".", 
    αβγ::AbstractVector{<:Real} = TEST_αβγ,
    isprimitive::Bool = true,
    flip_ksign::Bool = false
)

    sgnum === nothing && (sgnum = parse_sgnum(calcname))
    D === nothing     && (D = parse_dim(calcname))
    D < length(αβγ)   && (αβγ = αβγ[1:D])

    return _read_symdata(calcname, sgnum, Val(D), dir, αβγ, isprimitive, flip_ksign)
end

# function barrier cf. type-instability
function _read_symdata(
    calcname::AbstractString, 
    sgnum::Int, 
    Dᵛ::Val{D},
    dir::AbstractString, 
    αβγ::AbstractVector{<:Real},
    isprimitive::Bool, 
    flip_ksign::Bool
) where D
    
    # prepare default little groups of associated space group
    lgs⁰  = littlegroups(sgnum, Dᵛ)
    isprimitive && map!(g -> primitivize(g, #=modw=# false), values(lgs⁰)) # primitivize

    # read mpb dispersion data; use to guarantee frequency sorting at each k-point
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"), ',', Float64)::Matrix{Float64}
    kvs = collect(eachrow(@view dispersion_data[:,2:2+(D-1)]))
    Nk = length(kvs)
    freqs = dispersion_data[:,6:end]
    Nbands = size(freqs, 2)
    sortidxs = [sortperm(freqs_at_fixed_k) for freqs_at_fixed_k in eachrow(freqs)]

    # read mpb symmetry data, mostly as Strings (first column is Int)
    untyped_data = readdlm(joinpath(dir, calcname*"-symeigs.out"), ',', String; quotes=true)::Matrix{String}
    Nrows = size(untyped_data, 1)
    # process raw symmetry & operator data
    klabs    = Vector{String}(undef, Nk)                     # indexing: [kidx]
    lgs_ops  = Vector{Vector{SymOperation{D}}}(undef, Nk)    # indexing: [kidx][op]
    symeigsv = Vector{Vector{Vector{ComplexF64}}}(undef, Nk) # indexing: [kidx][band][op]
    rowidx = 1
    for kidx in 1:Nk
        kv = flip_ksign ? -kvs[kidx] : kvs[kidx] # `flip_ksign` is a hack to work around 
                                                 # phase convention issues; there be dragons
        klab = findfirst(lg->isapprox(position(lg)(αβγ), kv, atol=1e-6), lgs⁰)
        klab === nothing && error("could not find matching KVec for loaded kv = $kv")
        klabs[kidx]    = klab
        lgs_ops[kidx]  = Vector{SymOperation{D}}()
        symeigsv[kidx] = [Vector{ComplexF64}() for _ in 1:Nbands]
        while rowidx ≤ Nrows && parse(Int, untyped_data[rowidx, 1])::Int == kidx
            op = SymOperation{D}(strip(untyped_data[rowidx, 2], '"'))
            push!(lgs_ops[kidx], op)
            # symmetry eigenvalues at `op` (and `klab`) across all bands (frequency sorted)
            symeigs = parse.(ComplexF64, @view untyped_data[rowidx, 3:end][sortidxs[kidx]])
            push!.(symeigsv[kidx], symeigs) # broadcast over band indices
            rowidx += 1
        end
    end

    # build little groups
    lgs = map(zip(lgs_ops, klabs)) do (ops, klab)
        LittleGroup{D}(sgnum, position(lgs⁰[klab]), klab, ops)
    end

    return symeigsv, lgs
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return the irrep multiplicities for provided symmetry eigenvalues `symeigsv` with associated
little group irreps `lgirsv`, with symmetry eigenvalues aggregated (i.e. summed) over the 
band-indices in `bands`.

`symeigsv` and `lgirsv` must be `Vector`s referencing the same **k**-points, iterating in
shared order over the same little groups, and, subsequently, over the same little group
opeations. Additionally, the two must assume the same (primitive) setting.

Result is returned as a `Vector` over **k**-labels. If the eigenvalues at this **k**-point
are integer-expandable in the irreps in `lgirsv`, the `Vector`'s element is a vector of
integer multiplicities.
If not (e.g., due to `bands` "splitting" an irrep or overly tight `atol` tolerances),
the value is `nothing`.

## Keyword arguments
- `kwargs...`: optional keyword arguments (including `atol`) forwarded to 
  Crystalline.jl's `find_representation` (e.g., `atol` and `αβγ`).

If `bands` is unspecified, it assumes the maximal possible bandrange, as inferrable from the
size of `symeigsv`'s elements.

## Workflow
Symmetry eigenvalues `symeigsv` and little groups `lgs` can be loaded via 
[`read_symdata`](@ref); `lgirsv` can subsequently be generated via [`pick_lgirreps`](@ref).
To extract the multiplicities of the first 5 bands, we then call:
```jl
julia> bands = 1:5
julia> extract_multiplicities(symeigsv, lgirsv, bands)
```
"""
function extract_multiplicities( # ... main accessor
    symeigsv::Dict{String,<:Any},
    lgirsv::Dict{String,<:AbstractVector{LGIrrep{D}}},
    bands::AbstractVector{Int}=eachindex(first(values(symeigsv)));
    kwargs...
) where D
    if length(symeigsv) ≠ length(lgirsv)
        throw(DimensionError("mismatched lengths of `symeigsv` and `lgirsv`"))
    end
    multsv = Vector{Union{Nothing, Vector{Float64}}}(undef, length(symeigsv))
    for (kidx, (symeigs, lgirs)) in enumerate(zip(symeigsv, lgirsv))
        multsv[kidx] = find_representation(symeigs, lgirs, bands; kwargs...)
    end
    return multsv
end

"""
$(TYPEDSIGNATURES)

Return the irrep multiplicities for provided symmetry eigenvalues `symeigsv` with associated
little groups `lgs`, with symmetry eigenvalues aggregated (i.e. summed) over the 
band-indices in `bands`.

`symeigsv` and `lgs` must be `Vector`s referencing the same **k**-points and assuming the 
same operator sorting (and setting). `lgs` is used as an intermediary to determine the
relevant little group irreps via Crystalline.

Result is returned as a `Vector` over **k**-labels. If the eigenvalues at this **k**-point
are integer-expandable in the associated little group irreps, the `Vector` element is a
vector of multiplicities. If not (e.g., due to `bands` "splitting" an irrep or overly tight
`atol` tolerances), the value is `nothing`.

## Keyword arguments
- `timereversal` (default, `true`): whether the underlying MPB calculation has time-reversal
  symmetry; if `true`, physically real irreps (coreps) will be used. See Crystalline's
  `realify`.
- `isprimitive` (default, `true`): whether the underlying MPB calculation was performed in a
  primitive unit cell; if `true`, `lgs` is assumed to be provided in a primitive setting as
  well. Used as a consistency check in constructing associated little group irreps.
- `kwargs...`: optional keyword arguments (including `atol`) forwarded to Crystalline.jl's
  `find_representation` (e.g., `atol` and `αβγ`).

If `bands` is unspecified, it assumes the maximal possible bandrange, as inferrable from the
size of `symeigsv`'s elements.

## Workflow
Given a the "base" name of an MPB calculation, `calcname`, we can extract associated
multiplicities by:
```jl
julia> symeigsv, lgs = read_symdata(calcname, dir = "location/of/calcname/out/files")
julia> extract_multiplicities(symeigsv, lgs, 1:5)
```
where we consider the joint multiplicities of the first 5 bands. See also 
[`extract_all_multiplicities`](@ref).
"""
function extract_multiplicities( # ... main accessor
    symeigsv::AbstractVector,
    lgs::AbstractVector{LittleGroup{D}},
    bands=eachindex(first(values(symeigsv)));
    timereversal::Bool=true,
    isprimitive::Bool=true,
    kwargs...
) where D

    if length(symeigsv) ≠ length(lgs)
        error("`symeigsv` and `lgs` must have identical lengths")
    end
    
    lgirsv = pick_lgirreps(lgs; timereversal, isprimitive) # irreps for k-points in `lgs`

    return extract_multiplicities(symeigsv, lgirsv, bands; kwargs...)
end

function extract_multiplicities( # ... convenience accessor
    calcname::String,
    bands;
    timereversal::Bool=true,
    isprimitive::Bool=true,
    atol::Real = MULTIPLICITY_ATOL,
    αβγ::AbstractVector{<:Real} = TEST_αβγ,
    read_kwargs...
) 

    symeigsv, lgs = read_symdata(calcname; αβγ, isprimitive, read_kwargs...)
    return extract_multiplicities(symeigsv, lgs, bands; timereversal, isprimitive, atol, αβγ)
end

# ---------------------------------------------------------------------------------------- #

function pick_lgirreps(
    lgs::AbstractVector{LittleGroup{D}};
    timereversal::Bool=true,
    isprimitive::Bool=true
) where D
    sgnum = num(first(lgs))
    lgirsd = lgirreps(sgnum, Val(D))
    # retain the irreps corresponding to k-points in `lgs`
    lgirsv = [lgirsd[klabel(lg)] for lg in lgs]

    # incorporate timereveral & primitivize if relevant
    timereversal && map!(realify, lgirsv, lgirsv)
    if isprimitive
        cntr = centering(sgnum, D)
        if cntr ≠ 'P' || cntr ≠ 'p'
            for lgirs in lgirsv
                # all elements in `lgirs` point to the same group, so we should only change
                # _one_ such element - pick the first one (bit icky, yeah...)
                lg′ = primitivize(group(first(lgirs)), #=modw=# false)
                for (i, lgir) in enumerate(lgirs)
                    lgir′ = LGIrrep{D}(lgir.cdml, lg′, lgir.matrices, lgir.translations,
                                       lgir.reality, lgir.iscorep)
                    lgirs[i] = lgir′
                end
            end
        end
    end

    # ensure that `lgirsv` and `lgs` have identical operator sorting
    for (lgirs, lg) in zip(lgirsv, lgs)
        align_operators!(lgirs, lg) # mutates `lgirs` if operator sorting differs
    end

    return lgirsv
end

function align_operators!(lgirs::Vector{LGIrrep{D}}, lg::LittleGroup{D}) where D
    lg′ = group(first(lgirs))
    lg == lg′ && return lgirs
    length(lg) == length(lg′) || error("mismatched little groups")
    # find permutation to "realign" operator sorting between `lgs[kidx]` and `lgirsv[kidx]`
    perm = Vector{Int}(undef, length(lg))
    for (i, op) in enumerate(lg)
        idx = findfirst(==(op), lg′)
        idx === nothing && error("could not align operators")   
        perm[i] = idx
    end
    # apply permutation to `lgirs`
    for lgir in lgirs
        permute!.(lgir.matrices,     Ref(perm))
        permute!.(lgir.translations, Ref(perm))
    end
    permute!(first(lgirs).g.operations, perm) # NB: the `group` of each `lgirs` element is 
                        # actually the same struct, so we only permute it _once_
                        # (admittedly, a bit iffy to rely on...)
    @info "Permuted little group irreps to ensure aligned sorting" klab
    return lgirs
end
function align_operators!(lgirs::Collection{LGIrrep{D}}, lg::LittleGroup{D}) where D
    Collection(align_operators!(lgirs.vs, lg))
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return band-groupings and **k**-projected symmetry vectors at individual **k**-points
(`bandirsv`) and associated little group irreps (`lgirsv`).

To extract the associated **k**-point projected symmetry vectors of potentially separable
bands, see Crystalline.jl's `build_candidate_symmetryvectors`.

## Keyword arguments

- `latestarts` (default, `Dict("Γ" => D)`): allow user to specify that certain early bands
  be avoided for specific **k**-points. This is mainly useful to avoid considering the two
  singular Γ-point bands in 3D photonic crystals (which is the default behavior). For no
  skips, set to `nothing`. Defaults to the spatial dimension `D` at Γ (meaning, start at
  band `D` and skip the first `D-1` bands at Γ).
- `timereversal` & `isprimitive`: forwarded to [`read_symdata`](@ref) and 
  [`pick_lgirreps`](@ref).
- `atol` & `αβγ`: forwarded to Crystalline.jl's `find_multiplicities`.
- `read_kwargs`: additional keyword arguments forwarded to [`read_symdata`](@ref).

## Example

```jl
julia> calcname = "dim3-sg147-symeigs_6936-res32"
julia> bandirsv, lgirsv = 
    extract_all_multiplicities(calcname, timereversal=true, dir = "../../mpb-ctl/output/")
```
The result can be pretty-printed by e.g.:
```
julia> using Crystalline: label, symvec2string
julia> irlabs = [label.(lgirs) for lgirs in lgirsv]
julia> [[bands => symvec2string(n, irlabs[kidx]; braces=false) for (bands, n) in bandirs]
                                               for (kidx, bandirs) in enmerate(bandirsv)]
```
"""
function extract_all_multiplicities(
    calcname::String;
    timereversal::Bool = true,
    isprimitive::Bool = true,
    atol::Real = MULTIPLICITY_ATOL, 
    αβγ::AbstractVector{<:Real} = TEST_αβγ,
    latestarts::Dict{String,Int} = Dict("Γ" => parse_dim(calcname)),
    kwargs...
)
    symeigsv, lgs = read_symdata(calcname; αβγ, isprimitive, kwargs...)
    lgirsv = pick_lgirreps(lgs; timereversal, isprimitive)

    # determine how we can consistently group up irreps and bands across distinct k-points
    bandirsv = map(zip(lgirsv, symeigsv)) do (lgirs, symeigs)
        find_multiplicities(symeigs, lgirs; atol, αβγ, latestarts)
    end

    return bandirsv, lgirsv
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

Return the space group number of an MPB calculation `calcname` (in the conventions of
[`mpb_calcname`](@ref)).
"""
function parse_sgnum(calcname::AbstractString)
    idxs = findfirst("sg", calcname)
    idxs === nothing && error("could not infer space group number from calcname")
    # finds `sgnum` start & stop position in `calcname`, then pulls out `sgnum` as a string
    sgstart = idxs[end] + 1
    sgstop  = findnext(!isdigit, calcname, sgstart) - 1
    sgnum_str = calcname[sgstart:sgstop]
    
    return parse(Int, sgnum_str)
end

"""
$(TYPEDSIGNATURES)

Return the dimension of an MPB calculation `calcname` (in the conventions of
[`mpb_calcname`](@ref)).
"""
function parse_dim(calcname::AbstractString)
    idxs = findfirst("dim", calcname)
    idxs === nothing && error("could not infer dimension from calcname")

    Dstart = idxs[end] + 1
    D = parse(Int, calcname[Dstart]) # return the dimension as an integer
end
