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
B = matrix(brs)
F = smith(B)

# --- load and process data ---
bandirsd, lgirsd = extract_individual_multiplicities(
                        calcname;
                        timereversal,
                        dir = "../../mpb-ctl/output/sg147/")
length(lgirsd) ≠ length(klabels(sb)) && error("missing k-point data")

# extract the _potentially_ separable symmetry vectors `ns` and their band-ranges `bands`
bands, ns = extract_candidate_symmetryvectors(bandirsd, lgirsd, brs)
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
        is_transverse_bandstruct(n′, sb, lgirsd["Γ"], F)
    else                                                # regular bands
        isbandstruct(n′, F)
    end
    is_bs || continue

    # calculate topology of bands
    if first(band′) == 1
        topo = calc_topology_singular(n′, sb, lgirsd["Γ"], B)
    else
        topo = checkfragile ? calc_detailed_topology(n′, B, F) : calc_topology(n′, F)
    end

    println(band′, " ⇒ ", topo)
end