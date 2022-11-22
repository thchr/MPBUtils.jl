# returns a lazy "vector" of the (real) value of `flat` over the entire unit cell, using
# `nsamples` per dimension; avoids double counting at unit cell edges
function _lazy_fourier_eval_over_uc_as_vector(flat::AbstractFourierLattice{D}, nsamples::Integer) where D
    step = 1.0/nsamples
    samples = range(-0.5, 0.5-step, length=nsamples) # `-step` to avoid double counting ±0.5 (cf. periodicity)
    if D == 2
        itr = (real(flat(x,y)) for x in samples for y in samples)
    elseif D == 3
        itr = (real(flat(x,y,z)) for x in samples for y in samples for z in samples)
    end
end

"""
    $(TYPEDSIGNATURES)

Return the isovalue of `flat` such that the interior of the thusly defined levelset
isosurface encloses a fraction `filling` of the unit cell.

The keyword argument `nsamples` specifies the grid-resolution used in evaluating this
answer (via `quantile`).
"""
function filling2isoval(flat::AbstractFourierLattice{D}, filling::Real=0.5, nsamples::Integer=51) where D
    itr = _lazy_fourier_eval_over_uc_as_vector(flat, nsamples)
    return quantile(itr, filling)
end

"""
    $(TYPEDSIGNATURES)

Return the filling fraction of the interior of the isosurface defined by `flat` for an
isovalue `isoval`.

The keyword argument `nsamples` specifies the grid-resolution used in evaluating this
answer (using staircase integration).
"""
function isoval2filling(flat::AbstractFourierLattice{D}, isoval::Real, nsamples::Integer=51) where D
    itr = _lazy_fourier_eval_over_uc_as_vector(flat, nsamples)
    return count(<(isoval), itr)/nsamples^D
end
