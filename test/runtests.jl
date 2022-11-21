using GeneArrays
using Test

@testset "GeneArrays.jl" begin
    file = tempname()
    
    nrow, ncol = 200, 100
    a = rand(nrow, ncol)
    GeneArrays.writexy(file, a)
    b = GeneArrays.readxy(file)
    @test a == b

    a = a'a
    GeneArrays.writexy(file, a)
    b = GeneArrays.readxy(file)
    @test filesize(file) == ncol * (ncol + 1) * 4 + 23
    @test a == b

    rm(file)
end
