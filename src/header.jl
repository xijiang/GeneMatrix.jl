# © 2022, Xijiang Yu.  SPDX-License-Identifier: MIT

"""
This struct is designed for long term usage.  The idea magic header
was from plink bed file.  Matrix I/O is necessary as data nowadays can
be huge.  Disk I/O is not avoidable.
"""
struct xyheader
    x::Int8                     # X
    y::Int8                     # Y
    s::Int8                     # a space 0x20
    f::Int8                     # F L U S
    t::Int8                     # N T
    e::Int8                     # eltype
    r::Int8                     # '\n', reserved for future use
    m::Int64                    # nrow
    n::Int64                    # ncol
end

"""
    function mkhdr(nrow::Int64, ncol::Int64; mattp = 'F', trans = 'N', type = Int8)
Make a header for storage of a matrix in `XY` format.
"""
function mkhdr(nrow::Int64, ncol::Int64; 
               mattp = 'F',
               trans = 'N',
               type = Int8,
               )
    mattp ∈ "SUL" && nrow ≠ ncol && error("Matrix not square")
    code = findfirst(x -> x == type, vldtype)
    return xyheader('X', 'Y', ' ', mattp, trans, code, '\n', nrow, ncol)
end

"""
    function writehdr(oo::IOStream, hdr::xyheader)
Write an `xyheader` `hdr` to IOStream `oo`.
"""
function writehdr(oo::IOStream, hdr::xyheader)
    write(oo, [hdr.x, hdr.y, hdr.s, hdr.f, hdr.t, hdr.e, hdr.r])
    write(oo, [hdr.m, hdr.n])
end

"""
    function writehdr(file::AbstractString, hdr::xyheader)
Write a header `hdr` into `file`.  Matrix body can be appended after.
"""
function writehdr(file::AbstractString, hdr::xyheader)
    open(file, "w") do io
        writehdr(io, hdr)
    end
end

"""
    function readhdr(io::IOStream)
Read header from an IOStream in XY format.  Returns matrix storage
style, eltype and dimensions.
"""
function readhdr(io::IOStream)
    tmp = zeros(Int8, 7)
    read!(io, tmp)
    join(Char.(tmp[1:3])) == "XY " || error("Not in XY format")
    dim = zeros(Int64, 2)
    read!(io, dim)
    nrow, ncol = dim
    itp = tmp[6]
    itp > length(vldtype) && error("Type not defined")
    return Char(tmp[4]), Char(tmp[5]), vldtype[itp], nrow, ncol
end

"""
    function readhdr(file::AbstractString)
Fead header from a file of XY format.  Returns matrix storage style,
eltype and dimensions.
"""
function readhdr(file::AbstractString)
    isfile(file) || error("File $file not exist")
    open(file, "r") do io
        readhdr(io)
    end
end
