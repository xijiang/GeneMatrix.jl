using GeneArrays
using Test

@testset "I/O" begin
    file = tempname()
    f2   = tempname()

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

    a = rand(nrow, 2ncol)
    b = view(a, :, 1:ncol)
    c = view(a, :, ncol+1:2ncol)
    GeneArrays.writexy(file, b)
    GeneArrays.xyappend(file, c)
    d = GeneArrays.readxy(file)
    @test a == d

    GeneArrays.writexy(file, b)
    GeneArrays.writexy(f2, c)
    GeneArrays.xyappend(file, f2)
    d = GeneArrays.readxy(file)
    @test a == d

    rm(file)
    rm(f2)
end
