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
