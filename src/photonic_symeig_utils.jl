"""
$(TYPEDSIGNATURES)

Mutate the lowest `D-1` Γ-point entries of the provided vector of symmetry eigenvalues
`symeigsv` which correspond to the singular zero-frequency eigenstates. 
The mutated entries represent a valid, consistent choice of symmetry eigenvalues, following
the operator sorting of the little group (or little group irreps) in `lgs_or_lgirsv)`
associated with the Γ point. Returns the mutated `symeigsv`.

For `D = 2`, the polarization choice (transverse electric `:TE` or transverse magnetic
`:TM`) must be specified as the `polarization` argument.

## Details

The calculation of the singular zero-frequency symmetry eigenvalues of photonic crystals by
conventional numerical methods is usually unreliable since the photonic eigenvalue problem
is singular in the zero-frequency limit. To overcome this problem, the symmetry eigenvalues
of the `D-1` lowest bands at the Γ-point can be manually corrected, following the procedure
described in [^1].

**In 3D**, the two lowest-frequency symmetry eigenvalues ("2T") at the Γ-point are set such
that their sum is (``\\theta`` denoting the rotation angle of an operation ``g`` and
``det(g)`` its handedness):

``
x^{\\text{2T}}_{\\Gamma}(g) = \\det(g) (2\\cos(\\theta) + 1) - 1.
``

In practice, this sum is distributed arbitrarily over the two lowest bands at Γ.

**In 2D**, depending on the choice of `polarization` (either `:TM` or `:TE`), the
lowest-frequency Γ-point symmetry eigenvalue is set to:

``
x^{\\text{TM}}_{\\Gamma}(g) = 1, 
x^{\\text{TE}}_{\\Gamma}(g) = \\det(g).
``

[1]: Christensen, Po, Joannopoulos, & Soljacic, *Location and topology of the fundamental
     gap in photonic crystals*,
     [Physical Review X **12** 021066 (2022)](https://doi.org/10.1103/PhysRevX.12.021066).
"""
function fixup_gamma_symmetry!(
    symeigsv::AbstractVector{Vector{Vector{ComplexF64}}}, # ← mutated by call
    lgs_or_lgirsv::Union{AbstractVector{LittleGroup{D}}, AbstractVector{Collection{LGIrrep{D}}}},
    polarization::Union{Symbol, Nothing} = nothing
) where D
    Γ_idx = @something(
                findfirst(x->klabel(x)=="Γ", lgs_or_lgirsv),
                error("failed to find Γ-point"))
    lg_Γ = if lgs_or_lgirsv isa AbstractVector{LittleGroup{D}}
        lgs_or_lgirsv[Γ_idx]
    else # AbstractVector{Collection{LGIrrep{D}}}
        group(lgs_or_lgirsv[Γ_idx])
    end

    if D == 2
        if polarization == :TM
            symeigsv[Γ_idx][1] .= 1
        elseif polarization == :TE
            symeigsv[Γ_idx][1] .= det.(rotation.(lg_Γ))
        else
            throw(DomainError(polarization, "polarization must be :TE or :TM"))
        end
    elseif D == 3
        isnothing(polarization) || error("polarization argument must be `nothing` in 3D")
        ns = Crystalline.rotation_order.(lg_Γ) # rotation order (abs) and handedness (sign)
        x2T = sign.(ns) .* (2cospi.(2 ./ ns) .+ 1) .- 1
        splitting_fraction = 0.3721 # arbitrary fraction of 1 to avoid bands 1 and 2 
                                    # appearing splittable by accident
        symeigsv[Γ_idx][1] .= x2T .* splitting_fraction
        symeigsv[Γ_idx][2] .= x2T .* (1 - splitting_fraction)
    end
    return symeigsv
end

"""
    compute_symmetry_eigenvalues(ms, lgirsv[, polarization])

Given `ms`, an MPB `ModeSolver` `Py` object, and a vector of little group irrep tables,
`lgirsv`, sampled across distinct high-symmetry **k**-points, return the symmetry
eigenvalues across these little groups for the photonic bands returned by invoking `ms`.

The return value `symeigsv` is a triply nested `Vector`, whose elements
`symeigsv[kidx][n][i]` give the symmetry eigenvalue (or character) ⟨Eₙₖ|gᵢDₙₖ⟩ of the `n`th
band, for the `i`th little group operation gᵢ, in the `kidx`th little group (corresponding
to `group(lgirsv[kidx])`).

It is assumed that the `D-1` lowest frequency Γ point bands are zero-frequency solutions,
whose symmetry eigenvalues require manual correction. For 2D calculations, this correction
requires specification of the `polarization` argument (`:TE` or `:TM`). See also 
[`fixup_gamma_symmetry!`](@ref).
"""
function compute_symmetry_eigenvalues(
    ms::Py,
    lgirsv::AbstractVector{Collection{LGIrrep{D}}},
    polarization::Union{Nothing, Symbol} = nothing
) where D
    # compute and store symmetry eigenvalues ⟨Eₙₖ|gᵢDₙₖ⟩ in `symeigsv[kidx][n][i]`, indexing
    # into the `kidx´th little group (i.e., `lgirsv[kidx]`), the `n`th band index, and the 
    #`i`th operation (i.e., `group(lgirsv)[i]`)
    symeigsv = Vector{Vector{Vector{ComplexF64}}}(undef, length(lgirsv))
    N = pyconvert(Int, ms.num_bands)
    for (kidx, lgirs) in enumerate(lgirsv)
        lg = group(lgirs)
        kv = mp.Vector3(position(lg)()...)

        redirect_stdout(devnull) do # mute mpb's output
            ms.solve_kpoint(kv)
        end

        symeigsv[kidx] = [Vector{ComplexF64}(undef, length(lg)) for _ in 1:N]
        for (i, gᵢ) in enumerate(lg)
            # decompose gᵢ = {W|w}
            _W = rotation(gᵢ)
            W = if D == 3
                mp.Matrix(_W[:,1], _W[:,2], _W[:,3])
            elseif D == 2
                mp.Matrix(_W[:,1], _W[:,2], SVector(0.0, 0.0, 1.0))
            else # D == 1
                mp.Matrix(_W[:,1], SVector{3}(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
            end
            w = mp.Vector3(translation(gᵢ)...)
            symeigs = ms.compute_symmetries(W, w) # compute ⟨Eₙₖ|gᵢDₙₖ⟩ for all bands
            symeigs = pyconvert(Vector{ComplexF64}, symeigs) # convert Python to Julia types
            for (n, symeig) in enumerate(symeigs) # add to symeigs container
                symeigsv[kidx][n][i] = symeig
            end
        end
    end

    # --- fix singular photonic symmetry content at Γ, ω=0 ---
    fixup_gamma_symmetry!(symeigsv, lgirsv, polarization)

    return symeigsv
end