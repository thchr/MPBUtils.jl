using DelimitedFiles

# ---------------------------------------------------------------------------------------- #

struct BandSummary
    topology        :: TopologyKind
    band            :: UnitRange{Int}
    n               :: Vector{Int}
    brs             :: BandRepSet
    indicators      :: Vector{Int}
    indicator_group :: Vector{Int}
end

"""
$(TYPEDSIGNATURES)

Return the compatibility-allowed groupings of bands along with their topological properties
as a `Vector{BandSummary}`.

Keyword arguments `multiplicities_kwargs` are forwarded to
[`find_individual_multiplicities`](@ref).
"""
function analyze_symmetry_data(
            symeigsd::Dict{String, Vector{Vector{ComplexF64}}},
            lgirsd::Dict{String, Vector{LGIrrep{D}}},
            brs::BandRepSet;
            multiplicities_kwargs...) where D

    B = matrix(brs)
    F = smith(B)

    bandirsd = find_individual_multiplicities(symeigsd, lgirsd;
                                              multiplicities_kwargs...,
                                              latestarts = Dict{String, Int}())
    bands, ns = extract_candidate_symmetryvectors(bandirsd, lgirsd, brs;
                                                  latestarts = Dict{String, Int}())

    band_summaries = BandSummary[]
    idx = 1
    while idx ≤ length(ns)
        n_and_idx = _find_next_separable_band_grouping(ns, F, idx)
        if !isnothing(n_and_idx)
            # found a separable band grouping with sym vec `n′ = sum(ns[idx:idx′])`
            n′, idx′ = n_and_idx
            band′ = minimum(bands[idx]):maximum(bands[idx′])
            idx = idx′ + 1 # set to next "starting" index

            topo = calc_detailed_topology(n′, B; allow_nonphysical=true)
            band_summary = BandSummary(topo, band′, n′, brs,
                                       indicators(n′, F; allow_nonphysical=true)...)
            
            push!(band_summaries, band_summary)
        else
            break # could not find any more separable bands in `ns`; stop iteration
        end
    end
    return band_summaries
end

function _find_next_separable_band_grouping(
            ns::AbstractVector{<:AbstractVector{<:Integer}}, F::Smith, idx::Integer=1)

    idx > length(ns) && return nothing
    n′ = ns[idx]
    isbandstruct(n′, F; allow_nonphysical=true) && return n′, idx
    while (idx += 1) ≤ length(ns)
        n′ += ns[idx]
        isbandstruct(n′, F; allow_nonphysical=true) && return n′, idx
    end
    return nothing
end

function fix_gamma_irrep(symeigsd::Dict{String, Vector{Vector{ComplexF64}}}, lgd::Dict{String, LittleGroup{D}},
    polarization::Union{AbstractString, Nothing}=nothing) where D

    lg_gamma = lgd["Γ"]
    fix_gamma_irrep(symeigsd, lg_gamma, polarization)

end

function fix_gamma_irrep(symeigsd::Dict{String, Vector{Vector{ComplexF64}}}, lg_gamma::LittleGroup{D},
    polarization::Union{AbstractString, Nothing}=nothing) where D
    symeigsd_fixed = symeigsd
    signs = [det(rotation(g)) for g in lg_gamma]
    if D == 2
        if polarization == "tm"
            symeigsd_fixed["Γ"][1] = ones(length(signs))
        elseif polarization == "te"
            symeigsd_fixed["Γ"][1] = signs
        end
    elseif D ==3
        x2T = signs .* (2cospi.(2 ./ Crystalline.rotation_order.(lg_gamma)) .+ 1) .- 1
        symeigsd_fixed[Γ][1] = x2T*0.65
        symeigsd_fixed[Γ][2] = x2T*0.35
    end
    return symeigsd_fixed
end


Base.summary(io::IO, bs::BandSummary) = print(io, length(bs.band), "-band BandSummary:")
function Base.show(io::IO, ::MIME"text/plain", bs::BandSummary)
    summary(io, bs)
    println(io)
    println(io, " bands:      ", bs.band)
    print(io,   " n:          ", )
    Crystalline.prettyprint_symmetryvector(io, bs.n, irreplabels(bs.brs), braces=false)
    println(io)
    print(io, " topology:   ", lowercase(string(bs.topology)))
    if bs.topology == NONTRIVIAL
        print(io, "\n indicators: ")
        νs = bs.indicators
        length(νs) > 1 && print(io, "(")
        join(io, bs.indicators, ",")
        length(νs) > 1 && print(io, ")")
        printstyled(io, " ∈ ", classification(bs.indicator_group), color=:light_black)
    end
end
function Base.show(io::IO, bs::BandSummary) # compact print
    print(io, length(bs.band), "-band (", lowercase(string(bs.topology)), "): ")
    Crystalline.prettyprint_symmetryvector(io, bs.n, irreplabels(bs.brs))
end


"""
$(TYPEDSIGNATURES)

Return a vector of band ranges `bands` and associated symmetry vectors `ns` in the
irrep-sorting of `brs`, given a dictionary of symmetry data `bandirsd` (see
[`extract_individual_multiplicities`](@ref)) and a dictionary of correspoinding irrep data
`lgirsd`.

The returned symmetry vectors are _potentially_ separable, in the sense that they have
integer symmetry data and a consistent connectivity across **k**-points. Aside from this,
no account of compatibility relations are included. To test whether a returned symmetry
vector in fact fulfills compatibility relations, use SymmetryBases.jl's `is_bandstruct`.

## Keyword arguments
- `permd`: a dictionary of permutation vectors (indexed by **k**-vector labels), specifying
  a mapping between the irrep sortings in `lgirsd` and `brs` (see
  [`MPUtils.find_permutation`](@ref)). If set manually, `brs` should not be supplied (and
  will be ignored).

## Example
Supposing a set of symmetry data `bandirsd` has been extracted by
[`extract_individual_multiplicities`](@ref), 
"""
function extract_candidate_symmetryvectors(
            bandirsd::Dict{String, Vector{Pair{UnitRange{Int}, Vector{Int}}}},
            lgirsd::Dict{String, <:Vector{<:LGIrrep{D}}},
            brs=nothing;
            permd::Dict{String, Vector{Int}}=_default_permutation(lgirsd, brs),
            latestarts::Dict{String, Int}=Dict("Γ" => D)
            ) where D
    
    klabs = keys(bandirsd)
    length(lgirsd) ≠ length(klabs) && error("missing k-point data")
    
    bands, nds = collect_separable(bandirsd, lgirsd; latestarts)
    μs = length.(bands)
    isempty(bands) && error("found no isolable band candidates")
    # construct symmetry vectors, accounting for sorting mismatch specified by `permd`
    Nirs = sum(length, values(permd))
    ns = [Vector{Int}(undef, Nirs+1) for _ in 1:length(bands)]
    for (b, (nd, μ)) in enumerate(zip(nds, μs))
        for (klab, nᵏ) in nd
            permᵏ = permd[klab]
            ns[b][permᵏ] .= nᵏ
        end
        ns[b][end] = μ
    end
    return bands, ns
end

"""
$(TYPEDSIGNATURES)

Return a dictionary of permutation vectors, providing a mapping between the irrep sortings
and positions in `brs` relative to those in `lgirsd`.

## Example
```jl
julia> using Crystalline, MPBUtils

julia> brs = bandreps(210, 3);    # plane group 2

julia> lgirsd = lgirreps(210, 3);

julia> permd = find_permutation(lgirsd, brs) # indices of irreps in `lgirsd[klab]` in `brs`
Dict{String, Vector{Int64}} with 4 entries:
  "Y" => [1, 2]
  "B" => [3, 4]
  "A" => [5, 6]
  "Γ" => [7, 8]

julia> all([label.(lgirsd[klab]) == irreplabels(brs)[idxs] for (klab, idxs) in permd])
true
```
"""
function find_permutation(lgirsd::Dict{String, <:Vector{<:LGIrrep}}, brs::BandRepSet)
    return find_permutation(lgirsd, klabels(brs), irreplabels(brs))
end
function find_permutation(
            lgirsd::Dict{String, <:Vector{<:LGIrrep}},
            klabs::Vector{String},
            irlabs::Vector{String})
    # find permutation vectors between irreps in `lgirsd` and those in `irlabs` and `klabs`;
    # used to ensure alignment between sorting of symmetry vectors and band representations
    # since the alignment in `lgirsd` and e.g. `bandreps(...)` may differ
    permd = Dict(klab => Vector{Int}(undef, length(lgirsd[klab])) for klab ∈ klabs)
    for klab in klabs
    lgirs = lgirsd[klab]
        for (i, lgir) in enumerate(lgirs)
            irlab = label(lgir)
            j = findfirst(==(irlab), irlabs)
            j === nothing && error("Could not find irrep label $irlab")
            permd[klab][i] = j
        end
    end
    return permd
end

function _default_permutation(
            lgirsd::Dict{String, <:Vector{<:LGIrrep}},
            brs::Union{Nothing, BandRepSet})
    brs === nothing && error("must supply either `brs` or `permd`")
    return find_permutation(lgirsd, brs)
end

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
function read_symdata(calcname::AbstractString; 
            sgnum::Union{Int, Nothing}=nothing,
            D::Union{Int, Nothing}=nothing,
            dir::AbstractString=".", 
            αβγ::AbstractVector{<:Real}=TEST_αβγ,
            isprimitive::Bool=true,
            flip_ksign::Bool=false)

    sgnum === nothing && (sgnum = parse_sgnum(calcname))
    D === nothing     && (D = parse_dim(calcname))
    D < length(αβγ)   && (αβγ = αβγ[1:D])

    # prepare default little groups of associated space group
    lgs⁰  = littlegroups(sgnum, Val(D))
    isprimitive && map!(g -> primitivize(g, #=modw=# false), values(lgs⁰)) # primitivize

    # read mpb dispersion data; use to guarantee frequency sorting at each k-point
    dispersion_data = readdlm(joinpath(dir, calcname*"-dispersion.out"), ',')
    kvs = collect(eachrow(@view dispersion_data[:,2:2+(D-1)]))
    Nk     = length(kvs)
    freqs = dispersion_data[:,6:end]
    Nbands = size(freqs, 2)
    sortidxs = [sortperm(freqs_at_fixed_k) for freqs_at_fixed_k in eachrow(freqs)]
    # read mpb symmetry data, mostly as Strings (first column is Int)
    untyped_data = readdlm(joinpath(dir, calcname*"-symeigs.out"), ',', quotes=true)
    Nrows = size(untyped_data, 1)
    # process raw symmetry & operator data
    lgopsd   = Dict{String, Vector{SymOperation{D}}}()    # indexing: [klabel][op]
    symeigsd = Dict{String, Vector{Vector{ComplexF64}}}() # indexing: [klabel][band][op]
    rowidx = 1
    for kidx in 1:Nk
        kv = flip_ksign ? -kvs[kidx] : kvs[kidx] # `flip_ksign` is a hack to work around 
                                                 # phase convention issues; there be dragons
        klab = findfirst(lg->isapprox(position(lg)(αβγ), kv, atol=1e-6), lgs⁰)
        klab === nothing && error("could not find matching KVec for loaded kv = $kv")
        lgopsd[klab]   = Vector{SymOperation{D}}()
        symeigsd[klab] = [Vector{ComplexF64}() for _ in 1:Nbands]
        while rowidx ≤ Nrows && untyped_data[rowidx, 1] == kidx
            op = SymOperation{D}(strip(untyped_data[rowidx, 2], '"'))
            push!(lgopsd[klab], op)
            # symmetry eigenvalues at `op` (and `klab`) across all bands (frequency sorted)
            symeigs = parse.(ComplexF64, @view untyped_data[rowidx, 3:end][sortidxs[kidx]])
            push!.(symeigsd[klab], symeigs) # broadcast over band indices
            rowidx += 1
        end
    end
    # build little groups
    lgd = Dict(klab => LittleGroup{D}(sgnum, position(lgs⁰[klab]), klab, ops) for (klab, ops) in lgopsd)

    return symeigsd, lgd
end

# ---------------------------------------------------------------------------------------- #


function symeigs2irreps(symeigs::AbstractVector, lgirs::Vector{LGIrrep{D}}, bands;
            αβγ::AbstractVector{<:Real}=TEST_αβγ, kwargs...) where D
                                                                        # ... root accessor
    D < length(αβγ) && (αβγ = αβγ[1:D])
    
    # return multiplicities over provided irreps (and across `bands`)
    symeigs_bands = sum(@view symeigs[bands])
    return find_representation(symeigs_bands, lgirs, αβγ, Float64; kwargs...)
end

"""
$(TYPEDSIGNATURES) --> Dict{String, Union{Nothing, Vector{Float64}}}

Return the irrep multiplicities for provided symmetry eigenvalues `symeigsd` with associated
little group irreps `lgirsd`, with symmetry eigenvalues aggregated (i.e. summed) over the 
band-indices in `bands`.

`symeigsd` and `lgirsd` must be `Dict`s referencing the same **k**-points and assuming the 
same operator sorting (and setting).

Result is returned as a `Dict` over **k**-labels. If the eigenvalues at this **k**-point are
expandable in the irreps in `lgirsd`, the `Dict` value is a vector of multiplicities.
If not (e.g., due to `bands` "splitting" an irrep or overly tight `atol` tolerances),
the value is `nothing`.

## Keyword arguments
- `kwargs...`: optional keyword arguments (including `atol`) forwarded to 
  [`symeigs2irreps`](@ref) (e.g., `atol` and `αβγ`).

If `bands` is unspecified, it assumes the maximal possible bandrange, as inferrable from the
size of `symeigsd`'s elements.

## Workflow
Symmetry eigenvalues `symeigsd` and little groups `lgd` can be loaded via 
[`read_symdata`](@ref); `lgirsd` can subsequently be generated via [`pick__lgirreps`](@ref).
To extract the multiplicities of the first 5 bands, we then call:
```jl
julia> bands = 1:5
julia> extract_multiplicities(symeigsd, lgirsd, bands)
```
"""
function extract_multiplicities(symeigsd::Dict{String,<:Any}, 
            lgirsd::Dict{String,Vector{LGIrrep{D}}},
            bands=eachindex(first(values(symeigsd)));
            kwargs...) where D                                  # ... main accessor

    nd = Dict{String, Union{Nothing, Vector{Float64}}}()
    for klab in keys(symeigsd)
        nd[klab] = symeigs2irreps(symeigsd[klab], lgirsd[klab], bands; kwargs...)
    end
    return nd
end

"""
$(TYPEDSIGNATURES) --> Dict{String, Union{Nothing, Vector{Float64}}}

Return the irrep multiplicities for provided symmetry eigenvalues `symeigsd` with associated
little groups `lgd`, with symmetry eigenvalues aggregated (i.e. summed) over the 
band-indices in `bands`.

`symeigsd` and `lgd` must be `Dict`s referencing the same **k**-points and assuming the 
same operator sorting (and setting). `lgd` is used as an intermediary to determine the
relevant little group irreps via Crystalline.

Result is returned as a `Dict` over **k**-labels. If the eigenvalues at this **k**-point are
expandable in the associated little group irreps, the `Dict` value is a vector of
multiplicities. If not (e.g., due to `bands` "splitting" an irrep or overly tight `atol`
tolerances), the value is `nothing`.

## Keyword arguments
- `timereversal` (default, `true`): whether the underlying MPB calculation has time-reversal
  symmetry; if `true`, physically real irreps (coreps) will be used. See Crystalline's
  `realify`.
- `isprimitive` (default, `true`): whether the underlying MPB calculation was performed in a
  primitive unit cell; if `true`, `lgd` is assumed to be provided in a primitive setting as
  well. Used as a consistency check in constructing associated little group irreps.
- `kwargs...`: optional keyword arguments (including `atol`) forwarded to 
  [`symeigs2irreps`](@ref) (e.g., `atol` and `αβγ`).

If `bands` is unspecified, it assumes the maximal possible bandrange, as inferrable from the
size of `symeigsd`'s elements.

## Workflow
Given a the "base" name of an MPB calculation, `calcname`, we can extract associated
multiplicities by:
```jl
julia> symeigsd, lgd = read_symdata(calcname, dir = "location/of/calcname/out/files")
julia> extract_multiplicities(symeigsd, lgd, 1:5)
```
where we consider the joint multiplicities of the first 5 bands. See also 
[`extract_individual_multiplicities`](@ref).
"""
function extract_multiplicities(symeigsd::Dict{String, <:Any}, 
            lgd::Dict{String,LittleGroup{D}}, 
            bands=eachindex(first(values(symeigsd)));
            timereversal::Bool=true, isprimitive::Bool=true,
            kwargs...) where D                                  # ... main accessor

    if length(symeigsd) ≠ length(lgd) || !(keys(symeigsd) ⊆ keys(lgd))
        error("`symeigsd` and `lgd` must have identical keys")
    end
    
    lgirsd = pick_lgirreps(lgd; timereversal, isprimitive) # irreps for k-points in `lgd`

    return extract_multiplicities(symeigsd, lgirsd, bands; kwargs...)
end

function extract_multiplicities(calcname::String, bands;
            timereversal::Bool=true, isprimitive::Bool=true,
            atol::Real=DEFAULT_ATOL, αβγ::AbstractVector{<:Real}=TEST_αβγ,
            read_kwargs...)                                     # ... convenience accessor
    
    symeigsd, lgd = read_symdata(calcname; αβγ, isprimitive, read_kwargs...)
    return extract_multiplicities(symeigsd, lgd, bands; timereversal, isprimitive, atol, αβγ)
end

# ---------------------------------------------------------------------------------------- #

function pick_lgirreps(lgd::Dict{String, LittleGroup{D}};
            timereversal::Bool=true, isprimitive::Bool=true) where D
    sgnum = num(first(values(lgd)))
    lgirsd = lgirreps(sgnum, Val(D))
    # filter out k-points not in `lgd`
    filter!(((klab,_),) -> klab ∈ keys(lgd), lgirsd)

    # incorporate timereveral & primitivize if relevant
    if timereversal
        map!(realify, values(lgirsd))
    end
    if isprimitive
        cntr = centering(sgnum, D)
        if cntr ≠ 'P' || cntr ≠ 'p'
            for (klab, lgirs) in lgirsd
                # all elements in `lgirs` point to the same group, so we should only change
                # _one_ such element - pick the first one (bit icky, yeah...)
                lg = group(first(lgirs))
                lg.operations .= primitivize.(lg.operations, cntr, #=modw=# false)
            end
        end
    end

    # ensure that `lgirsd` and `lgd` have identical operator sorting
    for (klab, lg) in lgd
        align_operators!(lgirsd[klab], lg) # mutates `lgirsd` if operator sorting differs
    end

    return lgirsd
end

function align_operators!(lgirs::Vector{LGIrrep{D}}, lg::LittleGroup{D}) where D
    lg′ = group(first(lgirs))
    lg == lg′ && return lgirs
    length(lg) == length(lg′) || error("mismatched little groups")
    # find permutation to "realign" operator sorting between `lgd[klab]` and `lgirsd[klab]`
    perm = Vector{Int}(undef, length(lg))
    for op in lg
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

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES) --> bandirsd, lgirsd

Return band-groupings and **k**-projected symmetry vectors at individual **k**-points
(`bandirsd`) and associated little group irreps (`lgirsd`).

To extract the associated **k**-point projected symmetry vectors of potentially separable
bands, see [`collect_separable`](@ref) and [`merge_to_symvectors`](@ref).

## Keyword arguments

- `latestarts` (default, `Dict("Γ" => D)`): allow user to specify that certain early bands
  be avoided for specific **k**-points. This is mainly useful to avoid considering the two
  singular Γ-point bands in 3D photonic crystals (which is the default behavior). To skip 
  nothing, set to `Dict{String,Int}()`. Defaults to the spatial dimension `D` at Γ
  (meaning, start at band `D` and skip the first `D-1` bands at Γ).
- `timereversal` & `isprimitive`: forwarded to [`read_symdata`](@ref) and 
  [`extract_multiplicities`](@ref).
- `atol` & `αβγ`: forwarded to [`read_symdata`](@ref) and [`find_representation`](@ref).
- `read_kwargs`: additional keyword arguments forwarded to [`read_symdata`](@ref).

## Example

```jl
julia> calcname = "dim3-sg147-symeigs_6936-res32"
julia> bandirsd, lgirsd = extract_individual_multiplicities(calcname,
                        timereversal=true, dir = "../../mpb-ctl/output/", atol=1e-3)
```
The result can be pretty-printed by e.g.:
```
julia> using Crystalline: label, symvec2string
julia> irlabs = Dict(klab => label.(lgirs) for (klab, lgirs) in lgirsd)
julia> Dict(klab => [bands => symvec2string(n, irlabs[klab]; braces=false)
                     for (bands, n) in bandirs]         for (klab, bandirs) in bandirsd)
```
"""
function extract_individual_multiplicities(calcname::String;
            timereversal::Bool=true, isprimitive::Bool=true,
            atol::Real=DEFAULT_ATOL, αβγ::AbstractVector{<:Real}=TEST_αβγ,
            latestarts::Dict{String,Int}=Dict("Γ" => parse_dim(calcname)),
            kwargs...)

    symeigsd, lgd = read_symdata(calcname; αβγ, isprimitive, kwargs...)
    lgirsd = pick_lgirreps(lgd; timereversal, isprimitive)

    # determine how we can consistently group up irreps and bands across distinct k-points
    bandirsd = find_individual_multiplicities(symeigsd, lgirsd; atol, αβγ, latestarts)

    return bandirsd, lgirsd
end

function find_individual_multiplicities(symeigsd::Dict{String,<:AbstractVector},
            lgirsd::Dict{String,Vector{LGIrrep{D}}};
            atol::Real=1e-2,
            αβγ::AbstractVector{<:Real}=TEST_αβγ,
            latestarts::Dict{String,Int}=Dict("Γ" => D),
            maxresnorm::Real=1e-3) where D
    
    Nbands = length(first(values(symeigsd)))

    bandirsd = Dict(klab => Vector{Pair{UnitRange{Int}, Vector{Int}}}() for klab in keys(lgirsd))
    for (klab, lgirs) in lgirsd
        symeigs = symeigsd[klab]
        start = stop = get(latestarts, klab, 1)
        while stop ≤ Nbands
            bands = start:stop
            n = symeigs2irreps(symeigs, lgirs, bands; atol, αβγ, maxresnorm)
            if n !== nothing
                idxs = findall(nᵢ -> nᵢ > atol && abs(round(nᵢ) - nᵢ) < atol, n)
                if length(idxs) ≥ 1
                    # `bands` makes up a valid whole-multiple irrep combination at `klab`
                    # we emphasize that this can be a _combination_ i.e. multiple whole
                    # irreps; this can arise in cases of near-degeneracies, where numerical
                    # "mode-mixing" effectively spreads two irreps over two bands bands
                    # (even if they _should_ be each in one band); in principle, this can
                    # always be circumvented by increasing the resolution, but this is not
                    # a practical solution in general
                    push!(bandirsd[klab], bands => round.(Int, n))
                    start = stop + 1 # prepare for next band grouping
                end
            end
            stop += 1
        end
    end

    return bandirsd
end

"""
$(TYPEDSIGNATURES)

Return the "projected" symmetry-vectors for potentially separable bands using `bandirsd`
(see [`extract_individual_multiplicities`](@ref)) and `lgirsd` (see [`read_symdata`](@ref) 
and [`pick_lgirreps`](@ref)).
"""
function collect_separable(bandirsd::Dict{String, Vector{Pair{UnitRange{Int}, Vector{Int}}}},
                lgirsd::Dict{String, Vector{LGIrrep{D}}};
                latestarts::Dict{String, Int}=Dict("Γ" => D)) where D

    # Check for empty values in bandirsd and return gracefully. Otherwise the mapreduce
    # line will throw an error.
   if any(isempty, values(bandirsd))
      # return empty `collectibles_bands, collectibles_symvecs`
      return Vector{UnitRange{Int}}(), Vector{Dict{String, Vector{Int}}}()
   end

    Nbands = mapreduce(last∘first∘last, min, values(bandirsd)) # smallest "last" band-index
    include = Dict(klab => Int[] for klab in keys(bandirsd))
    collectibles_bands   = Vector{UnitRange{Int}}()
    collectibles_include = Vector{typeof(include)}()
    start = stop = 1
    while stop ≤ Nbands
        stable = true
        for (klab, bandirs) in bandirsd
            idxs = include[klab]
            latestart = get(latestarts, klab, nothing)
            for (i, (bands, _)) in enumerate(bandirs)
                minband, maxband = extrema(bands)
                minband < start && continue

                if latestart !== nothing && stop < latestart && minband ≤ latestart
                    break # allow `bands` to "not count" if `latestarts` indicates skips
                end
                i ∉ idxs && push!(idxs, i)
                
                if maxband > stop
                    stop = maxband
                    stable = false
                    break
                elseif maxband == stop
                    break
                end
            end
        end

        if stable
            push!(collectibles_bands,   start:stop)
            push!(collectibles_include, include)
            start = (stop += 1)
            include = Dict(klab => Int[] for klab in keys(bandirsd))
        end
    end
    
    # create projected symmetry vectors for each, using knowledge of number of irreps in
    # each little group
    Nirs = Dict(klab => length(lgirsd[klab]) for klab in keys(bandirsd))
    collectibles_symvecs = Vector{Dict{String, Vector{Int}}}()
    for include in collectibles_include
        symvecs = Dict(klab => zeros(Int, Nir) for (klab, Nir) in Nirs)
        for (klab, idxs) in include
            for i in idxs
                symvecs[klab] += last(bandirsd[klab][i])
            end
        end
        push!(collectibles_symvecs, symvecs)
    end

    return collectibles_bands, collectibles_symvecs
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
