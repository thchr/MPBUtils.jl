"""
    BandSummary{D} <: Crystalline.AbstractSymmetryVector{D}

A summary of the band symmetry and topology of a set of bands. `BandSummary`s of consecutive
bands can be added, corresponding to stacking of bands.
"""
struct BandSummary{D} <: Crystalline.AbstractSymmetryVector{D}
    topology        :: SymmetryBases.TopologyKind
    bands           :: UnitRange{Int}
    n               :: SymmetryVector{D}
    brs             :: Collection{NewBandRep{D}}
    indicators      :: Vector{Int}
    indicator_group :: Vector{Int}
end

# ---------------------------------------------------------------------------------------- #
# AbstractSymmetryVector interface

Crystalline.SymmetryVector(bs::BandSummary) = bs.n
Base.size(bs::BandSummary) = size(bs.bands)

# ---------------------------------------------------------------------------------------- #
# arithmetic

function Base.:+(bs1::BandSummary{D}, bs2::BandSummary{D}) where D
    # check bands have almost certainly the same bandreps (but only check cheap things)
    num(bs1.brs) == num(bs2.brs) || error("bands had band representations from different space groups")
    first(bs1.brs).timereversal == first(bs2.brs).timereversal || error("bands had band representations with different time-reversal symmetry")
    length(bs1.brs) == length(bs2.brs) || error("bands had band representations with different cardinality")

    # consecutiveness check for `band`
    if last(bs1.bands) + 1 == first(bs2.bands)
        band = first(bs1.bands):last(bs2.bands)
    elseif last(bs2.bands) + 1 == first(bs1.bands)
        band = first(bs2.bands):last(bs1.bands)
    else
        throw(DomainError((bs1.bands, bs2.bands), "bands of bs1 and bs2 must be consecutive"))
    end
    n = bs1.n + bs2.n
    indicators = mod.(bs1.indicators .+ bs2.indicators, bs1.indicator_group)
    topology = !iszero(indicators) ? NONTRIVIAL :
                    SymmetryBases.calc_detailed_topology(n, bs1.brs; allow_negative=true)

    return BandSummary{D}(topology, band, n, bs1.brs, indicators, bs1.indicator_group)
end

function Base.:(==)(bs1::BandSummary{D}, bs2::BandSummary{D}) where D
    return bs1.n == bs2.n && bs1.bands == bs2.bands
end

# ---------------------------------------------------------------------------------------- #
# printing

function Base.summary(io::IO, bs::BandSummary{D}) where D
    print(io, length(bs.bands), "-band BandSummary{", D, "}")
end

function Base.show(io::IO, ::MIME"text/plain", bs::BandSummary)
    summary(io, bs)
    println(io, ":")
    println(io, " bands:      ", bs.bands)
    println(io, " n:          ", bs.n)
    print(io,   " topology:   ")
    t = bs.topology
    printstyled(io, lowercase(string(t));
                color = t==NONTRIVIAL ? :green : t==FRAGILE ? :yellow : :light_red)
    if bs.topology == NONTRIVIAL
        print(io, "\n indicators: ")
        νs = bs.indicators
        length(νs) > 1 && print(io, "(")
        join(io, bs.indicators, ",")
        length(νs) > 1 && print(io, ")")
        printstyled(io, " ∈ ", indicator_group_as_string(bs.indicator_group),
            color=:light_black)
    end
end

function Base.show(io::IO, bs::BandSummary) # compact print
    t = bs.topology
    print(io, bs.n, " ")
    printstyled(io, lowercase(string(t)); 
                color = t==NONTRIVIAL ? :green : t==FRAGILE ? :yellow : :light_red)
end

# ---------------------------------------------------------------------------------------- #

"""
$(TYPEDSIGNATURES)

An extended variant of Crystalline's `collect_compatible`, which additionally performs a
symmetry-indicator/TQC-based topological analysis for each band, returning a vector of
`BandSummaries{D}` with the associated information.
"""
function collect_compatible_detailed(
    symeigsv::AbstractVector{Vector{Vector{ComplexF64}}},
    brs::Collection{NewBandRep{D}},
    B::AbstractMatrix{<:Integer} = stack(brs),
    F::Smith{<:Integer} = smith(B);
    kws...
) where D
    ns = collect_compatible(symeigsv, brs, F; kws...) # get symmetry vectors
    b = 1
    summaries = Vector{BandSummary{D}}(undef, length(ns))
    inds_group = indicator_group(brs)
    for (i, n) in enumerate(ns)
        b′ = b+occupation(n)-1
        bands = b:b′
        inds = symmetry_indicators(n, F; allow_negative=true)
        topo = if iszero(inds)
            SymmetryBases.calc_detailed_topology(n, B; allow_negative=true)
        else
            NONTRIVIAL
        end
        summaries[i] = BandSummary{D}(topo, bands, n, brs, inds, inds_group)
        b = b′+1
    end
    return summaries
end