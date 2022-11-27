module GeneMatrix

using LinearAlgebra, DataFrames

abstract type TwoBM end
# 1:14. New type must be appened in the end to keep the package consistant.
const vldtype = (Bool, Int8, Int16, Int32, Int64, Int128, UInt8, UInt16,
                 UInt32, UInt64, UInt128, Float16, Float32, Float64)

# Write your package code here.
include("header.jl")
include("io.jl")
include("nrm.jl")               # numerical relationship matrix
include("grm.jl")               # genomic relationship matrix
include("graph.jl")        # store huge genomic data with graph theory
include("pedigree.jl")          # pedigree related operations

export readxy, writexy, readbed, writebed
end
