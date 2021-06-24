using Crystalline
using Crystalline: formatirreplabel
using SymmetryBases
using PhotonicBandConnectivity
using MPBUtils

# --- setup info ---
sgnum  = 147
has_tr = false
id     = 7716
checkfragile = true
calcname = "dim3-sg147-breaktr-detfix-g1.0_symeigs_"*string(id)*"-res32"

# --- hilbert/ebr bases ---
sb, brs = compatibility_basis(sgnum, 3; timereversal=has_tr)
B = matrix(brs)
F = smith(B)

# --- load and process data ---
bandirsd, lgirsd = extract_individual_multiplicities(
                        calcname,
                        timereversal=has_tr,
                        dir = "../../mpb-ctl/output/sg147/",
                        atol=2e-2)
length(lgirsd) ≠ length(sb.klabs) && error("missing k-point data")

bands, nds = collect_separable(bandirsd, lgirsd)
μs = length.(bands)

isempty(bands) && error("   ... found no isolable band candidates ...")

# --- build symmetry vectors ---
# find the permutation between sorting in `sb` and `lgirsd`
permd = Dict(klab => Vector{Int}(undef, length(lgirsd[klab])) for klab ∈ sb.klabs)
for klab in sb.klabs
    lgirs = lgirsd[klab]
    for (i, lgir) in enumerate(lgirs)
        irlab = formatirreplabel(label(lgir))
        j = findfirst(==(irlab), sb.irlabs)
        j === nothing && error("Could not find irrep label $irlab")

        permd[klab][i] = j
    end
end

# use permutation to construct symmetry vectors, in `sb`'s sorting
ns = [Vector{Int}(undef, length(first(sb))) for _ in 1:length(bands)]
for (b, (nd, μ)) in enumerate(zip(nds, μs))
    for (klab, nᵏ) in nd
        permᵏ = permd[klab]
        ns[b][permᵏ] .= nᵏ
    end
    ns[b][end] = μ
end

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
        topo = topology_from_2T1L_xor_1L(n′, sb, lgirsd["Γ"], B)
    else
        topo = checkfragile ? calc_detailed_topology(n′, B, F) : calc_topology(n′, F)
    end

    println(band′, " ⇒ ", topo)
end