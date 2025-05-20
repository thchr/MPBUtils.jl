using Crystalline
using MPBUtils
write_dir = (@__DIR__)*"" # NB: make sure whatever directory this points to actually exists

function init_lattice(sgnum, Dᵛ::Val{D}, cntr::Char;
                maxGs=ntuple(_->2, Val(D)), expon::Real=1.25) where D
    # level-set surface
    flat = normscale!(modulate(levelsetlattice(sgnum, Dᵛ, maxGs)), expon)
    deleteat!(flat.orbits, 1), deleteat!(flat.orbitcoefs, 1) # trivial constant G=0 term

    return  primitivize(flat, cntr) # return a primitive lattice
end

# --- choices ---
write_to_file = true
D = 3       # dimension
sgnum = 68  # space group number
tr = true   # whether we have time-reversal or not
cntr = centering(sgnum, D) # centering symbol (e.g, 'F' for face-centered, etc.)
res = 32    # resolution used in MPB
nbands = 12 # number of bands requested from MPB

# --- find out which little groups we need to assess bandreps/symmetry vector ---
brs = primitivize(calc_bandreps(sgnum, D, timereversal=tr))
lgs = group.(irreps(brs))

# --- generate a bunch of mpb-input files ---
ids = 1:10 # create 10 different lattices and "input-files"
runtype = "all" # we don't have TE/TM in 3D; do all polarizations
for id in ids
    # direct basis
    Rs  = directbasis(sgnum, Val(D)) # generate a random basis consistent with `sgnum`
    pRs = primitivize(Rs, cntr)

    # generate a isosurface
    pflat = init_lattice(sgnum, Val(D), cntr)

    # generate some dielectric values
    filling = rand(0.25:0.05:0.75)
    epsin = rand((7.0, 9.0, 11.0, 13.0, 15.0))
    epsout = 1.0

    # write to .sh file for mpb
    filename = mpb_calcname(D, sgnum, id, res, runtype)
    filepath = joinpath(write_dir, filename*".sh")
    write_to_file || continue
    open(filepath, "w") do io
        prepare_mpbcalc!(io, sgnum, pflat, pRs, filling, epsin, epsout, runtype;
                             res=res, lgs=lgs, id=id, nbands=nbands)
        # NB: if you want an mpb dispersion calculation, provide a keyword argument `kvs`
        #     instead of `lgs` in the above. You can also supply a filename for `kvs` if
        #     you want.
    end
    id == 1 && println("MPB setup files written to:")
    println("  ∙  ", filepath, ".sh")
end