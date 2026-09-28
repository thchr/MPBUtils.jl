using Crystalline, MPBUtils, Brillouin, GLMakie, ProgressMeter, PythonCall
using LinearAlgebra: norm
mp = pyimport("meep")
mpb = pyimport("meep.mpb")

# --- mpb: geometry & solver initialization ---
D = 3
sgnum = 225 # space group number
Rs′ = directbasis(sgnum, Val(3))
Rs = primitivize(Rs′, centering(sgnum))

m = mp.Medium(epsilon=1)
f = 1
geometry = [mp.Sphere(center=mp.Vector3(0,0,0), radius=1/sqrt(8)*f, material=m)]
lattice = mp.Lattice(basis_size=norm.(Rs), # take units relative to conventional unit cell
                     basis1 = Rs[1], basis2 = Rs[2], basis3 = Rs[3])
ms = mpb.ModeSolver(
    num_bands        = 20,
    k_points         = [],
    geometry         = pylist(geometry),
    geometry_lattice = lattice,
    resolution       = 16,
    tolerance        = 1e-6,
    default_material = mp.Medium(epsilon=13),
)
ms.init_params(p = mp.NO_PARITY, reset_fields = true)

# --- band representations, littlegroups, & irreps ---
brs = primitivize(bandreps(sgnum, Val(D))) # elementary band representations
lgirsv = irreps(brs)                       # associated little groups & small irreps

# --- compute band symmetry data ---
symeigsv = compute_symmetry_eigenvalues(ms, lgirsv)

# --- analyze connectivity and topology of symmetry data ---
summaries = collect_compatible_detailed(symeigsv, brs)

# --- band structure ---
kp = irrfbz_path(sgnum, Rs′);
kvs = interpolate(kp, 150);

freqs = Matrix{Float64}(undef, length(kvs), pyconvert(Int, ms.num_bands));
@showprogress 0.1 for (i, kv) in enumerate(kvs)
    redirect_stdout(devnull) do
        ms.solve_kpoint(mp.Vector3(kv...))
    end
    freqs[i,:] = sort!(pyconvert(Vector{Float64}, ms.get_freqs()))
end

plot(kvs, freqs; annotations = collect_irrep_annotations(symeigsv, lgirsv))