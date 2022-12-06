using GeneMatrix
using Test

@testset "I/O" begin
    file = tempname()
    f2   = tempname()

    nrow, ncol = 200, 100
    a = rand(nrow, ncol)
    GeneMatrix.writexy(file, a)
    b = GeneMatrix.readxy(file)
    @test a == b

    a = a'a
    GeneMatrix.writexy(file, a)
    b = GeneMatrix.readxy(file)
    @test filesize(file) == ncol * (ncol + 1) * 4 + 24
    @test a == b

    a = rand(nrow, 2ncol)
    b = view(a, :, 1:ncol)
    c = view(a, :, ncol+1:2ncol)
    GeneMatrix.writexy(file, b)
    GeneMatrix.xyappend(file, c)
    d = GeneMatrix.readxy(file)
    @test a == d

    GeneMatrix.writexy(file, b)
    GeneMatrix.writexy(f2, c)
    GeneMatrix.xyappend(file, f2)
    d = GeneMatrix.readxy(file)
    @test a == d
    
    rm(file)
    rm(f2)
end
