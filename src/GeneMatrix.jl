module GeneMatrix

using LinearAlgebra

# 1:14. New type must be appened in the end to keep the package consistant.
const vldtype = (Bool, Int8, Int16, Int32, Int64, Int128, UInt8, UInt16,
                 UInt32, UInt64, UInt128, Float16, Float32, Float64)

# Write your package code here.
include("header.jl")
include("mats.jl")
include("io.jl")

end
