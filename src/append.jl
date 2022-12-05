"""
    xyappend(file::AbstractString, mat::Matrix; trans='N')

Append genotypes in `mat` to `file`.
"""
function xyappend(file::AbstractString, mat::Matrix; trans='N')
    isfile(file) || error("File $file doesn't exist")
    open(file, "r+") do io
        mt, _, et, nrow, ncol = readhdr(io)
        mt ∈ "LUS" && error("Can't merge triangles")
        et ≠ eltype(mat) && error("Matrices have different eltypes")

        seekend(io)
        if trans == 'N'
            ncl2 = size(mat, 2)
            nrow ≠ size(mat, 1) && error("Matrices don't match")
            write(io, mat)
        else
            ncl2 = size(mat, 1)
            nrow ≠ size(mat, 2) && error("Matrices don't match")
            write(io, mat')
        end

        seek(io, 15)
        write(io, [ncol + ncl2])
    end
end

"""
    xyappend(fa::AbstractString, fb::AbstractString)

Append genotypes in `mb` to file `ma`.  This is only valid if the
matrix type is `F`, and two matrices have same number of ID, or same
number of loci.
"""
function xyappend(fa::AbstractString, fb::AbstractString)
    (isfile(fa) && isfile(fb)) || error("File $fa or $fb doesn't exist")

    mtx, ntx, etx, nrwx, nclx = readhdr(fa)
    mty, nty, ety, nrwy, ncly = readhdr(fb)
    (mtx == 'F' && mty == 'F') || error("Can't merge triangles")
    (etx == ety) || error("Not of t`he same eltype matrices")
    trans = (ntx == nty) ? 'N' : 'T'

    mat = Mmap.mmap(fb, Matrix{ety}, (nrwy, ncly), 24)
    xyappend(fa, mat, trans = trans)
end
