# ---------------------------------------------------------------------------------------- #
# ⚠ LEGACY EXAMPLE — kept for archival purposes; not a suggested usage pattern
#
# This is the analysis half of the old two-step workflow: its counterpart
# `legacy-example-scm-setup-script.jl` writes MPB input files, MPB is then run separately,
# and this script reads the resulting `-dispersion.out` and `-symeigs.out` files back in.
#
# For new work, drive MPB directly from Julia instead and use `compute_symmetry_eigenvalues`
# together with `collect_compatible_detailed`; see `inverse-opal.jl` and the README. That
# path needs no intermediate files, and it keeps the little groups, operator sorting, and
# phase conventions consistent by construction.
#
# The zero-frequency (singular) band analysis below additionally relies on
# PhotonicBandConnectivity.jl, whose 2T+1L treatment requires an explicit choice of
# longitudinal symmetry vector `nᴸ`; see that package for what that choice means.
# ---------------------------------------------------------------------------------------- #

using Crystalline
using SymmetryBases
using PhotonicBandConnectivity
using MPBUtils

# --- setup info ---
D = 3
sgnum = 147
timereversal = false
id = 7716
checkfragile = true
calcname = "dim3-sg147-breaktr-detfix-g1.0_symeigs_"*string(id)*"-res32"

# --- hilbert/ebr bases ---
sb, brs = compatibility_basis(sgnum, D; timereversal)
B = stack(brs)
F = smith(B)

# --- load and process data ---
bandirsv, lgirsv = extract_all_multiplicities(
                        calcname;
                        timereversal,
                        dir = "../../mpb-ctl/output/sg147/")
length(lgirsv) ≠ length(klabels(sb)) && error("missing k-point data")

idx_Γ = something(findfirst(lgirs->klabel(lgirs)=="Γ", lgirsv))
lgirs_Γ = lgirsv[idx_Γ]

# the singular ω=0 analysis needs an explicit longitudinal symmetry vector `nᴸ`; take the
# one picked by PhotonicBandConnectivity (`nothing` if the space group requires no 1L pick)
_, _, sb¹ᴸ, idx¹ᴸ = minimal_expansion_of_zero_freq_bands(sgnum; timereversal)
nᴸ = isnothing(idx¹ᴸ) ? nothing : sb¹ᴸ[idx¹ᴸ]

# extract the _potentially_ separable symmetry vectors `ns` and their band-ranges `bands`
ns = Crystalline.build_candidate_symmetryvectors(bandirsv, lgirsv;
                                                 latestarts = Dict("Γ" => D))
stops = cumsum(occupation.(ns))
bands = UnitRange.([1; stops[1:end-1] .+ 1], stops)

# --- find which band combinations are in {BS} and then test their topology ---
band′ = 0:0
n′ = similar(first(ns))
is_bs = true # at every new iteration, `is_bs` effectively means `was_prev_iter_bs`
for (band, n) in zip(bands, ns)
    global band′, n′, is_bs

    band′ = is_bs ? band : (minimum(band′):maximum(band))
    n′    = is_bs ? n : n′ + n

    # test if `n′` is in {BS} or not
    is_bs = if first(band′) == 1                        # singular bands
        is_transverse_bandstruct(Vector(n′), sb, lgirs_Γ, F)
    else                                                # regular bands
        iscompatible(n′, F)
    end
    is_bs || continue

    # calculate topology of bands
    if first(band′) == 1
        isnothing(nᴸ) && (println(band′, " ⇒ (no 1L pick; skipping)"); continue)
        topo = calc_topology_singular(Vector(n′), nᴸ, F)
    else
        topo = checkfragile ? calc_detailed_topology(n′, B, F) : calc_topology(n′, F)
    end

    println(band′, " ⇒ ", topo)
end
