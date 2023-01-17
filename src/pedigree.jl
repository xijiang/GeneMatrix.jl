"""
    codeped(id, pa, ma)
Code these 3 string vectors into integers starts from 1.
A DataFrame is returned which has 3 columns,
`sire::Int`, `dam::Int`. and `name::AbstractString`.
String "0" and number 0 is deemed as unknown.
Row numbers of the DataFrame are the coded ID number.
The `Name` column contains the original ID names.

The reason to have 3 String vectors here is to make this as general
as possible.  One needs to prepare them for this function,
which is easy.
"""
function codeped(id, pa, ma)
    # QC
    length(id) == length(pa) == length(ma) || 
        error("Lengths of id, pa, ma differ")
    length(unique(id)) == length(id) || error("Repeated ID")
    op = setdiff(setdiff(pa, id), ["0"])
    om = setdiff(setdiff(ma, id), ["0"])
    isempty(intersect(op, om)) || error("ID is/are both sire and dam")

    df = DataFrame(sire = Int32[], dam = Int32[], name = AbstractString[])
    dc = Dict{AbstractString, Int}()
    dc["0"] = uid = 0

    for x in union(op, om)
        uid += 1
        dc[x] = uid
        push!(df, (0, 0, x))
    end

    remain = ref = length(id)
    while ref > 0
        for i in eachindex(id)
            haskey(dc, id[i]) && continue
            haskey(dc, pa[i]) || continue
            haskey(dc, ma[i]) || continue
            uid += 1
            dc[id[i]] = uid
            push!(df, (dc[pa[i]], dc[ma[i]], id[i]))
            remain -= 1
        end
        ref > remain || error("Looped pedigree")
        ref = remain
    end
    df
end
