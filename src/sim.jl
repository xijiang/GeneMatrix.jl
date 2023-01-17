"""
    quickBase(nlc::Int, nid::Int; maf = .05, bp = .75)
A quick way to simulate genotypes of `nLoci` by `nID`.
Allele frequencies are sampled from a `Beta(.75, .75)`,
conditioned on `maf`.
It has a U-shaped distribution by default.
"""
function quickBase(nlc::Int, nid::Int; maf = .05, bp = .75)
    0 < maf < 0.5 || error("Maf $maf not in (0, 0.5)")
    gt = zeros(Int8, nlc, nid)
    for iic in 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(Beta(bp, bp))
        end
        rand!(Binomial(2, p), view(gt, iic, :))
    end
    gt
end

"""
    function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
Normalize QTL effect, such that the TBV variance is within `1.0 ± ϵ`.
"""
function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
    nqtl, nid = size(Q)
    bv = zeros(nid)
    matmul!(bv, Q', efct)
    m, s = mean(bv), std(bv)
    while abs(m) > ϵ || abs(s - 1) > ϵ
        efct .-= m/nqtl
        efct ./= s
        matmul!(bv, Q', efct)
        m, s = mean(bv), std(bv)
    end
end

"""
    function simQTL(gt::Matrix{Int8}, nqtl...; d = Laplace(), norm = true, ϵ = 1e-5)
## Description
Given genotype `gt` of `nLoc × nID`, this function sample `nqtl` loci as QTL.
QTL effects are distributed as Laplace by default.  After this sampling, their 
signs are randomly flipped.
This also normalize the true breeding values to have approximate mean 0, 
and variance 1.

The function returns an array of named tupples.

## Example distribution
- `Laplace()`, ≡ `Laplace(0, 1)`
- `Gamma()`, ≡ `Gamma(1, 1)`
- `Normal()`, ≡ `Normal(0, 1)`.
"""
function simQTL(gt::AbstractArray, nqtl...; d = Laplace(), norm = true, ϵ = 1e-5)
    nlc, nid = size(gt)
    qtl = []
    for n in nqtl
        loci = sort(randperm(nlc)[1:n])
        efct = rand(d, n) .* rand([-1, 1], n)
        Q = copy(view(gt, loci, :))
        norm && norm_qtl(Q, efct, ϵ)
        p = vec(mean(Q, dims = 2) ./ 2)
        v = 2 .* p .* (1 .- p) .* efct .^ 2
        push!(qtl, (locus = loci,
                    effect = efct,
                    ernk = sortperm(abs.(efct), rev = true),
                    vrnk = sortperm(v, rev = true)
                    ))
    end
    qtl
end

#=
"""
    function simQTL(fgt::String, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
File layer of function `simQTL`.
"""
function simQTL(fgt::String, nqtl...; d = Laplace(), ϵ = 1e-5, norm = true)
    nlc, nid = Fio.readdim(fgt)
    gt = Mmap.mmap(fgt, Matrix{Int8}, (nlc, nid), 24)
    qtl = simQTL(gt, nqtl, d = d, norm = norm, ϵ = ϵ)
end
=#

"""
    function simPtQTL(gt, nqtl; d = MvNormal(zeros(2), I(2)))
Simulation `nqtl` pleiotropic QTL for two traits, with genotype `gt` of `nLoci × nID`.
"""
function simPtQTL(gt, nqtl; d = MvNormal(zeros(2), I(2)))
    nlc, nid = size(gt)
    loci = sort(randperm(nlc)[1:nqtl])
    efct = rand(d, nqtl)
    return (locus = loci, e1 = efct[1, :], e2 = efct[2, :])
end

"""
    simMap(len, nlc...; cM = 1_000_000)

Simulate a linkage map DataFrame with chromosome number of ``UInt8`` and 
base pair positions of ``UInt32``.
This can be used for crossover simulations.
"""
function simMap(len, nlc...; cM = 1_000_000)
    lmp = DataFrame(chr = UInt8[], pos = UInt32[])
    ich = UInt8(1)
    for n in nlc
        pos = sort(rand(1:len, n))
        append!(lmp, DataFrame(chr=repeat([ich], length(pos)), pos = pos))
        ich += 1
    end
    lmp
end
