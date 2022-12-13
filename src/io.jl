# © 2022, Xijiang Yu.  SPDX-License-Identifier: MIT

"""
    function readxy(file)
Read a matrix `mat` from `file` according to the specifications in its
header.
"""
function readxy(file)
    isfile(file) || error("File $file doesn't exist")
    filesize(file) < 23 &&
        error("File $file too small to be an 'XY' matrix")
    mat = nothing               # to return if success
    open(file, "r") do io
        mt, id, et, nrow, ncol = readhdr(io)
        if et == Bool           # eltype
            mat = BitArray(undef, nrow, ncol)
            read!(io, mat)
        else
            mat = zeros(et, nrow, ncol)
            if mt == 'L' || mt == 'S' # matrix type
                for i in 1:ncol
                    read!(io, view(mat, i:nrow, i))
                end
                if mt == 'S'
                    for i in 1:ncol
                        copyto!(view(mat, i, i+1:ncol), view(mat, i+1:nrow, i)')
                    end
                end
            elseif mt == 'U'
                for i in 1:ncol
                    read!(io, view(mat, 1:i, i))
                end
            else
                read!(io, mat)
            end
                    
        end
    end    
    return mat
end

"""
    function writexy(file, mat; mattp = 'F', trans = 'N')
Write a matrix `mat` into `file`, with specified saving type.  Note a
symmetric matrix is written of its lower triangle.
"""
function writexy(file, mat; mattp = 'F', trans = 'N')
    mattp ∈ "FLUS" || error("Matrix type $mattp not defined")
    trans ∈ "NT" || error("ID locus type $trans not defined")
    et = eltype(mat)
    issymmetric(mat) && (mattp = 'S') # force write half of a symm. mat.
    nrow, ncol = size(mat)
    
    hdr = mkhdr(nrow, ncol, mattp = mattp, trans = trans, type = et)
    open(file, "w") do io
        writehdr(io, hdr)
        if mattp == 'L' || mattp == 'S'
            for i in 1:ncol
                write(io, mat[i:nrow, i])
            end
        elseif mattp == 'U'
            for i in 1:ncol
                write(io, mat[1:i, i])
            end
        else
            write(io, mat)
        end
    end
end

"""
    function readbed(bed, nid)
Read genotypes in file `bed` of `nid` samples and return a 2-bit matrix.
For reference, the two-bit values and their meaning:

| Value | Meaning |
| :---: | :------ |
| 00 | Homozygous for first allele in .bim file |
| 01 | Missing genotype |
| 10 | Heterozygous |
| 11 | Homozygous for second allele in .bim file |

"""
function readbed(bed, nid)
    isfile(bed) || error("File $bed doesn't exist")
    fs = filesize(bed)
    bs = Int(ceil(nid/4))
    @show bs
    nlc = (fs - 3) ÷ bs
    nlc * bs + 3 == fs || error("File $bed has wrong size")
    mat = nothing
    open(bed, "r") do io
        magic = Int8[0x6c, 0x1b, 0x01]
        tst = copy(magic)
        read!(io, magic)
        magic == tst || error("File $bed has no magic header")
        mat = zeros(Int8, bs, nlc)
        read!(io, mat)
    end
    mat
end

"""
    function writebed(bed, mat)
Write Int8 matrix `mat` into plink bed format.
"""
function writebed(bed, mat)
    magic = Int8[0x6c, 0x1b, 0x01]
    eltype(mat) == Int8 || error("Can't write non Int8 matrix to plink bed")
    open(bed, "w") do io
        write(io, magic)
        write(io, mat)
    end
end

function readmacs(dir)
end
