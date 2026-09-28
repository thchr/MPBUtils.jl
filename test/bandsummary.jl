using Crystalline
using MPBUtils
using LinearAlgebra: tr
using Test

# build the exact symmetry eigenvalues `symeigsv[kidx][band][op]` of a symmetry vector `n`,
# spreading each irrep's characters evenly over the bands that the irrep occupies
function synthetic_symeigs(n::Crystalline.AbstractSymmetryVector)
    map(zip(irreps(n), multiplicities(n))) do (lgirs, mults)
        symeigs = Vector{Vector{ComplexF64}}()
        for (lgir, m) in zip(lgirs, mults)
            iszero(m) && continue
            χ = tr.(lgir()) # characters of `lgir` over `group(lgir)`
            d = Crystalline.irdim(lgir)
            append!(symeigs, fill(χ ./ d, m*d))
        end
        symeigs
    end
end

@testset "BandSummary" begin
for (sgnum, D) in ((10, 2), (13, 2), (81, 3), (230, 3))
    brs = primitivize(bandreps(sgnum, Val(D)))
    @test brs isa Collection{BandRep{D, LGIrrep{D}, SiteIrrep{D}}}

    # a single EBR must come back as a single, trivial band grouping of equal content
    br = brs[1]
    summaries = collect_compatible_detailed(synthetic_symeigs(br), brs)
    bs = only(summaries)
    @test bs isa BandSummary{D}
    @test SymmetryVector(bs) == SymmetryVector(br)
    @test bs.bands == 1:occupation(br)
    @test bs.topology == TRIVIAL
    @test iszero(bs.indicators)

    # a `BandSummary` must behave as the symmetry vector it wraps
    @test length(bs) == length(SymmetryVector(br))
    @test collect(bs) == collect(SymmetryVector(br))
    @test occupation(bs) == occupation(br)

    # stacking must preserve the total symmetry content
    n = SymmetryVector(brs[1]) + SymmetryVector(brs[min(2, length(brs))])
    summaries = collect_compatible_detailed(synthetic_symeigs(n), brs)
    @test sum(SymmetryVector, summaries) == n
    @test sum(occupation, summaries) == occupation(n)
    if length(summaries) ≥ 2
        bs′ = summaries[1] + summaries[2]
        @test bs′ isa BandSummary{D}
        @test SymmetryVector(bs′) == SymmetryVector(summaries[1]) +
                                     SymmetryVector(summaries[2])
        @test bs′.bands == first(summaries[1].bands):last(summaries[2].bands)
    end
end
end # @testset
