"""
    fixup_bloch_phases!(symeigsv, lgs_or_lgirsv, αβγ = nothing)

Convert the symmetry eigenvalues in `symeigsv`, as computed by MPB, from MPB's Bloch phase
convention to that of Crystalline.jl's little group irreps, and return the mutated
`symeigsv`.

`symeigsv[kidx][n][i]` must give the symmetry eigenvalue of band `n` under the `i`th
operation of the little group associated with `lgs_or_lgirsv[kidx]` (a `LittleGroup` or a
`Collection{LGIrrep}`), in the same operator sorting. `αβγ` sets the free parameters of any
non-special **k**-vectors.

## Details

MPB and Crystalline.jl associate opposite phases with the translation part of a symmetry
operation. For an operation ``g = \\{R|\\mathbf{w}\\}`` in the little group of
``\\mathbf{k}``, MPB returns the symmetry eigenvalue in the usual physical convention, i.e.,
that obtained by letting ``g`` act on Bloch states ``\\psi_{n\\mathbf{k}}(\\mathbf{r}) =
\\mathrm{e}^{\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{r}}u_{n\\mathbf{k}}(\\mathbf{r})`` as
``(g\\psi)(\\mathbf{r}) = \\psi(g^{-1}\\mathbf{r})``:

``
\\chi^{\\text{MPB}}_{n\\mathbf{k}}(g) = \\mathrm{e}^{-\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{w}} \\chi^{0}_{n\\mathbf{k}}(R),
``

with ``\\chi^{0}`` the part that depends only on the rotation ``R``. Crystalline.jl's
`LGIrrep`s instead follow the convention picked by ISOTROPY, the Bilbao Crystallographic
Server, and the Inui et al. textbook, namely, ``D^{\\mathbf{k}}(\\{R|\\mathbf{w}\\}) =
\\mathrm{e}^{+\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{w}}D^{\\mathbf{k}}(\\{R|\\mathbf{0}\\})``, i.e., they
expect

``
\\chi^{\\text{Crystalline}}_{n\\mathbf{k}}(g) = \\mathrm{e}^{+\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{w}} \\chi^{0}_{n\\mathbf{k}}(R),
``

corresponding to an assumption of Bloch waves propagating with a phase
``\\mathrm{e}^{-\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{r}}``.

The conversion is therefore a per-operation phase:

``
\\chi^{\\text{Crystalline}}_{n\\mathbf{k}}(g) = \\mathrm{e}^{+2\\mathrm{i}\\mathbf{k}\\cdot\\mathbf{w}} \\chi^{\\text{MPB}}_{n\\mathbf{k}}(g).
``

The correction is the identity unless ``2\\mathbf{k}\\cdot\\mathbf{w} \\notin \\mathbb{Z}``,
so it never does anything for symmorphic operations, nor at Γ. Even where the phase is
nontrivial it often affects no change on the inferred irrep multiplicities (e.g., because
the associated symmetry eigenvalues vanish).[^1]

[1]: But this is not always the case: at the P point of space groups 88 and 230, e.g., the
     uncorrected symmetry eigenvalues do not decompose into the little group irreps at all,
     and the analysis would silently return a wrong irrep assignment without this
     conversion.

The sign convention discrepancy is the same one documented in Crystalline.jl (issue #12).
"""
function fixup_bloch_phases!(
    symeigsv::AbstractVector{<:AbstractVector{<:AbstractVector{<:Number}}},
    lgs_or_lgirsv::Union{AbstractVector{LittleGroup{D}},
                         AbstractVector{Collection{LGIrrep{D}}}},
    αβγ::Union{AbstractVector{<:Real}, Nothing} = nothing
) where D
    length(symeigsv) == length(lgs_or_lgirsv) ||
        throw(DimensionMismatch("mismatched lengths of `symeigsv` and `lgs_or_lgirsv`"))
    for (symeigs, lg_or_lgirs) in zip(symeigsv, lgs_or_lgirsv)
        lg = lg_or_lgirs isa LittleGroup{D} ? lg_or_lgirs : group(lg_or_lgirs)
        kv = position(lg)(αβγ)
        # NB: Note that there is an additional factor of 2π here cf. the reduced coordinate
        #     representations of `kv` and `op`
        phases = [cispi(4*dot(kv, translation(op))) for op in lg]
        all(isone, phases) && continue # nothing to do at this k-point
        for symeigs_n in symeigs
            length(symeigs_n) == length(lg) ||
                error("`symeigsv` and its little group must agree on the operator count")
            symeigs_n .*= phases
        end
    end
    return symeigsv
end
