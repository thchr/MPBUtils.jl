const IntOrNothing = Union{Int, Nothing}

# ---------------------------------------------------------------------------------------- #
"""
$(TYPEDSIGNATURES)
"""
function load_symdata(calcname::AbstractString, 
                      sgnum::IntOrNothing=nothing, D::IntOrNothing=nothing; parentdir::AbstractString="./", 
                      αβγ::AbstractVector{<:Real}=Crystalline.TEST_αβγ,
                      flip_ksign::Bool=false)
    sgnum === nothing && (sgnum = parse_sgnum(calcname))
    D === nothing && (D = parse_dim(calcname))
    if D < length(αβγ)
        αβγ = αβγ[Base.OneTo(D)]
    end
    dispersion_data = readdlm(parentdir*calcname*"-dispersion.out", ',') 
    kvecs = KVec.(eachrow(@view dispersion_data[:,2:2+(D-1)])) 
    kidxs = eachindex(kvecs)
    Nk = length(kidxs)

    freqs = dispersion_data[:,6:end] 
    ordering_perms = [sortperm(freqs_at_fixed_k) for freqs_at_fixed_k in eachrow(freqs)]
    # read mpb data, mostly as Strings (first column is Int)
    untyped_data = readdlm(parentdir*calcname*"-symeigs.out", ',', quotes=true)
    Nrows = size(untyped_data, 1)
    # convert to typed format, assuming "Int, String, (ComplexF64 ...)" structure
    ops = [Vector{SymOperation{D}}() for _ in kidxs] # indexed first across KVecs, then across SymOperations
    symeigs =[Vector{Vector{ComplexF64}}() for _ in kidxs] # (--||--), then across bands
    kidx = rowidx = 1
    while kidx ≤ Nk
        while rowidx ≤ Nrows && untyped_data[rowidx, 1] == kidx
            op = SymOperation{D}(strip(untyped_data[rowidx, 2], '"'))
            push!(ops[kidx], op)
            # symmetry eigenvalues (frequency-sorted) at current `kidx` and current `op`
            _symeigs = parse.(ComplexF64, @view untyped_data[rowidx, 3:end][ordering_perms[kidx]])
            push!(symeigs[kidx], _symeigs)
            rowidx += 1
        end
        kidx += 1
    end
    # build little groups and reconstruct labels
    lgs⁰  = [primitivize(lg) for lg in values(get_littlegroups(sgnum, Val(D)))] # primitivizes & converts dict to vector...
    idxs⁰ = if !flip_ksign
        [findfirst(lg->isapprox(kvec(lg)(αβγ), constant(kv), atol=1e-6), lgs⁰) for kv in kvecs]
    else # for trickery, to change the interpreted sign of k
        [findfirst(lg->isapprox(kvec(lg)(αβγ), -constant(kv), atol=1e-6), lgs⁰) for kv in kvecs]
    end
    klabs = klabel.(getindex.(Ref(lgs⁰), idxs⁰))
    lgs = [LittleGroup{D}(sgnum, kvecs[kidx], klabs[kidx], ops[kidx]) for kidx in kidxs]
    return lgs, symeigs
end


function find_degeneracies(calcname::String; dir::String="./")
    dispersion_data = readdlm(dir*calcname*"-dispersion.out", ',') 
    D=2
    _, numbands = size(dispersion_data[:,6:end])
    possible_degeneracies = Vector{Vector{Vector{Integer}}}()
    for kvfreq in eachrow(dispersion_data[:,6:end])
        possible_degeneraciesk = Vector{Vector{Integer}}()
        rounded_freqs = round.(kvfreq, digits=3)
        unique_freqs = unique(rounded_freqs)
        for (index, uniquefreq) in enumerate(unique_freqs)
            atuniquefreq = Vector{Integer}()
            for i in 1:numbands
                rounded_freqs[i] ≈ uniquefreq ? push!(atuniquefreq, i) : nothing
            end
            length(atuniquefreq)>1 ? push!(possible_degeneraciesk, atuniquefreq) : nothing
        end
        push!(possible_degeneracies, possible_degeneraciesk)
    end
    return possible_degeneracies
end

"""
$(TYPEDSIGNATURES)

Return the symmetry data associated with `calcname` and the associated `LGIrrep`s for bands
summed over `bandidxs`.
If not supplied explicitly via `sgnum` and `D`, the method will attempt to infer the space group number
and dimension from `calcname`.

# Keyword arguments
- `timereversal`: whether to use real (`true`, default) or complex irreps (`false`).
- `isprimitive`: If the lattice in `calcname` is *not* in a primitive setting, this kwarg
  must be set to `false` (this information is needed as part of a check that ensures operator 
  equivalence across the little groups in `calcname` and those on from `get_lgirreps`).
- `atol`: absolute tolerance passed to `find_representation` ($(Crystalline.DEFAULT_ATOL),
  default)
- `αβγ`: value used for pinning the little groups of lines, planes, and volumes to a
  concrete k-vector. *Must* be identical to that used in setting up the calculation (e.g.,
  as in `write_lgs_to_mpb`, where it defaults to `Crystalline.TEST_αβγ` as it
  also does here).
  This setting is immaterial for "special" little groups, i.e. little groups without free
  k-vector parameters.

# Example
```julia-repl
using Crystalline: prettyprint_symmetryvector, formatirreplabel
using MPBUtils

calcname  = "dim3-sg68-symeigs_15-res32"
bandidxs  = 1:2
msvec, lgirsvec = symdata2representation(calcname, bandidxs)
msvec_int = [convert.(Int, round.(ms, digits=3)) for ms in msvec]
n         = vcat(msvec_int...) # symmetry vector
irlabs    = formatirreplabel.(label.(Iterators.flatten(lgirsvec)))

prettyprint_symmetryvector(stdout, n, irlabs)
```
Prints: `[Γ₂⁺+Γ₄⁺, T₁, Y₂⁺+Y₂⁻, Z₂, R₁, S₁]`.

**Note:** Γ irreps are not well-defined in this example (as they touch ``ω = 0``).
"""
function symdata2representation(calcname::AbstractString, bandidxs::AbstractVector=1:2,
            sgnum::IntOrNothing=nothing, D::IntOrNothing=nothing; parentdir::AbstractString="./", 
            timereversal::Bool=true, isprimitive::Bool=true, 
            atol::Float64=Crystalline.DEFAULT_ATOL,
            αβγ::AbstractVector{<:Real}=Crystalline.TEST_αβγ,
            lgidxs::Union{Nothing, AbstractVector{<:Integer}}=nothing,
            flip_ksign::Bool=false)
    sgnum === nothing && (sgnum = parse_sgnum(calcname))
    D === nothing && (D = parse_dim(calcname))
    if D < length(αβγ)
        αβγ = αβγ[Base.OneTo(D)]
    end
    lgs, symeigs = load_symdata(calcname, sgnum, D; parentdir=parentdir, αβγ=αβγ, flip_ksign=flip_ksign)
    if lgidxs !== nothing # possibly only look at a subset of all the little groups
        symeigs = symeigs[lgidxs]
        lgs     = lgs[lgidxs]
    end
    lgirsd = get_lgirreps(sgnum, D)
    lgirsvec = [lgirsd[klabel(lg)] for lg in lgs]
    timereversal && (lgirsvec .= realify.(lgirsvec))
    # check: operator sorting and values are consistent across lgs and lgirsvec
    lgs′ = group.(first.(lgirsvec))
    isprimitive && lgs′ .= primitivize.(lgs′, #=modw=#false)
    @assert lgs == lgs′
    # compute matching irreps at each k-point, i.e. the symmetry vector
    symvals = map(kidx->sum.(getindex.(symeigs[kidx], Ref(bandidxs))), eachindex(lgs))
    msvec   = map(kidx->find_representation(symvals[kidx], lgirsvec[kidx], αβγ, Float64; atol=atol), eachindex(lgs))
    return msvec, lgirsvec
end


function symdata2representation(lgs::Array{LittleGroup{D},1}, symeigs::Vector{<:Vector{<:Vector{<:Complex{Float64}}}}, bandidxs::AbstractVector=1:2,
    timereversal::Bool=true, isprimitive::Bool=true, 
    atol::Float64=Crystalline.DEFAULT_ATOL,
    αβγ::AbstractVector{<:Real}=Crystalline.TEST_αβγ,
    lgidxs::Union{Nothing, AbstractVector{<:Integer}}=nothing) where D

    sgnum = first(lgs).num 
    println("Dimension is ", D)
    if lgidxs !== nothing # possibly only look at a subset of all the little groups
        symeigs = symeigs[lgidxs]
        lgs     = lgs[lgidxs]
    end
    lgirsd = get_lgirreps(sgnum, D)
    lgirsvec = [lgirsd[klabel(lg)] for lg in lgs]
    timereversal && (lgirsvec .= realify.(lgirsvec))
    # check: operator sorting and values are consistent across lgs and lgirsvec
    lgs′ = group.(first.(lgirsvec))
    isprimitive && lgs′ .= primitivize.(lgs′, #=modw=#false)
    @assert lgs == lgs′
    # compute matching irreps at each k-point, i.e. the symmetry vector
    symvals = map(kidx->sum.(getindex.(symeigs[kidx], Ref(bandidxs))), eachindex(lgs))
    msvec   = map(kidx->find_representation(symvals[kidx], lgirsvec[kidx], αβγ, Float64; atol=atol), eachindex(lgs))
    return msvec, lgirsvec
end

function symdata2representation(lg::LittleGroup{D}, symeig::Vector{<:Vector{<:Complex{Float64}}}, bandidxs::AbstractVector=1:2,
    #=
    Only give one little group and one symeig to obtain the msvec corresponding to that particular kvector
    =#
    timereversal::Bool=true, isprimitive::Bool=true, 
    atol::Float64=Crystalline.DEFAULT_ATOL,
    αβγ::AbstractVector{<:Real}=Crystalline.TEST_αβγ) where D
    sgnum = lg.num
    lgirsd = get_lgirreps(sgnum, D)
    lgir = lgirsd[klabel(lg)] 
    timereversal && (lgir = realify(lgir))
    # check: operator sorting and values are consistent across lgs and lgirsvec
    lg′ = group(first(lgir))
    isprimitive && (lg′ = primitivize(lg′, #=modw=#false))
    @assert lg == lg′
    symvals = sum.(getindex.(symeig, Ref(bandidxs)))
    msvec   = find_representation(symvals, lgir, αβγ, Float64; atol=atol)
    return msvec, lgir
end

"""
$(TYPEDSIGNATURES)
Returns the irreps of a single band in the order of kvectors supplied in the MPB calculation. 
"""
function singlebandirreps(calcname::AbstractString, bandidx::Integer; parentdir::AbstractString="./")
    irrepsofband = Vector{Vector{String}}()
    msvec, lgirsvec = symdata2representation(calcname, bandidx:bandidx, parentdir=parentdir)
    irreplabels = [label.(lgirsatk) for lgirsatk in lgirsvec]
    for (index, (ms, lgirs, irreplabel)) in enumerate(zip(msvec, lgirsvec, irreplabels))
        try 
            msint = Int.(round.(ms, digits=3))
            @assert length(msint) == length(lgirs)
            @assert length(msint) == length(irreplabel)
            push!(irrepsofband, irreplabel[findall(x->!(x≈0), msint)])
        catch InexactError
            println("Likely found fractional index- checking for degeneracies")
            degeneratevec = find_degeneracies(calcname, dir=parentdir)[index]
            fractional_indices = findall(x-> !isinteger(round(x, digits=3)), ms)
            degeneratebands = degeneratevec[(findfirst(x->x>0, bandidx .∈ degeneratevec))]
            msint = round.(symdata2representation(calcname, degeneratebands)[1][index], digits=3)
            push!(irrepsofband, irreplabel[findall(x->!(x≈0), msint)])
            println(irreplabel[fractional_indices]...)
        end
    end
    return irrepsofband
end

"""
$(TYPEDSIGNATURES)
Return the spacegroup index of an MPB calculation by conventions set by mpb_calcname
"""
function parse_sgnum(calcname::AbstractString)
    sgstart = findfirst("sg", calcname)[end] + 1  #Finds sg in the string, selects the last index corresponding to sg, and then selects the index after that. 
    sgstop  = findnext(!isdigit, calcname, sgstart) - 1 
    #equivalent to findnext(x->isdigit(x)==false, calcname, sgstart )-1 
    sgnum   = parse(Int, calcname[sgstart:sgstop])  
    #We find where the digits start and stop and interpret the intervening characters as an integer
end

"""
$(TYPEDSIGNATURES)
Return the dimensionality of an MPB calculation by conventions set by mpb_calcname
"""
function parse_dim(calcname::AbstractString)
    Dstart = findfirst("dim", calcname)[end] + 1 #findfirst returns the range of characters that correspond to dim. +1 corresponds to the where the dimension is canonically set by mpb_calcname
    D = parse(Int, calcname[Dstart]) # return the dimension as an integer
end
