using Crystalline
using MPBUtils
using LinearAlgebra: dot, tr
using Test

# the inverse of `fixup_bloch_phases!`: take characters in Crystalline's convention to MPB's
function to_mpb_convention(χ, lg)
    kv = position(lg)()
    χ .* [cispi(-4*dot(kv, translation(op))) for op in lg]
end

@testset "Bloch phase conventions" begin
for sgnum in (88, 230, 81, 10, 13)
    D = sgnum ≤ 17 ? 2 : 3
    lgirsv = irreps(primitivize(bandreps(sgnum, Val(D))))
    for lgirs in lgirsv
        lg = group(lgirs)
        # characters of each irrep, in Crystalline's convention, one "band" per irrep
        χs = [tr.(lgir()) for lgir in lgirs]
        symeigsv = [[to_mpb_convention(χ, lg) for χ in χs]]

        fixup_bloch_phases!(symeigsv, [lg])
        @test only(symeigsv) ≈ χs

        # and each recovered character must decompose onto exactly its own irrep
        for (j, χ) in enumerate(only(symeigsv))
            @test find_representation(χ, lgirs) == [i == j for i in eachindex(lgirs)]
        end
    end
end

# the correction only ever touches operations with a fractional translation, and only where
# `2k⋅w ∉ ℤ`: it must be a no-op at Γ, and for every symmorphic group
for (sgnum, D) in ((230, 3), (88, 3), (10, 2))
    lgirsv = irreps(primitivize(bandreps(sgnum, Val(D))))
    for lgirs in lgirsv
        lg = group(lgirs)
        symeigsv = [[ComplexF64[i for i in eachindex(lg)]]]
        untouched = copy(only(only(symeigsv)))
        fixup_bloch_phases!(symeigsv, [lg])
        if klabel(lg) == "Γ" || sgnum == 10 # symmorphic ⇒ no fractional translations
            @test only(only(symeigsv)) == untouched
        end
    end
end

# the correction must undo `to_mpb_convention` exactly, and must do something somewhere
lgirsv = irreps(primitivize(bandreps(230, Val(3))))
original = [[ComplexF64[i + im*j for i in eachindex(group(lgirs))] for j in 1:3]
            for lgirs in lgirsv]
symeigsv = [[to_mpb_convention(χ, group(lgirs)) for χ in χs]
            for (χs, lgirs) in zip(original, lgirsv)]
@test symeigsv ≉ original # at P, where `2k⋅w ∉ ℤ`
@test fixup_bloch_phases!(symeigsv, lgirsv) ≈ original
end # @testset
