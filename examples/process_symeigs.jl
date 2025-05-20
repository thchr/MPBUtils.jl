using Crystalline
using SymmetryBases
using PhotonicBandConnectivity
using MPBUtils

# --- setup info ---
sgnum = 147
timereversal = false
id = 7716
checkfragile = true
calcname = "dim3-sg147-breaktr-detfix-g1.0_symeigs_"*string(id)*"-res32"

# --- hilbert/ebr bases ---
sb, brs = compatibility_basis(sgnum, 3; timereversal)
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

# extract the _potentially_ separable symmetry vectors `ns` and their band-ranges `bands`
bands, ns = Crystalline.build_candidate_symmetryvectors(bandirsv, lgirsv, brs)
μs = last.(ns) # associated connectivities

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
        is_transverse_bandstruct(n′, sb, lgirs_Γ, F)
    else                                                # regular bands
        iscompatible(n′, F)
    end
    is_bs || continue

    # calculate topology of bands
    if first(band′) == 1
        topo = calc_topology_singular(n′, sb, lgirs_Γ, B)
    else
        topo = checkfragile ? calc_detailed_topology(n′, B, F) : calc_topology(n′, F)
    end

    println(band′, " ⇒ ", topo)
end