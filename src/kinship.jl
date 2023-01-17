"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
The first 2 columns of ped must be `pa` and `ma`.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + .5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return .5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return .5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end

"""
    function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int, Int}, Float64})
Recursive kinship calculation with kinship of ID pair `(i, j)`
stored in dictionary `dic`.  The first 2 columns of ped must be `pa` and `ma`.
The memory usage may be bigger than Meuwissen and Luo 1992, or Quaas 1995.
The speed is however standable.
The recursive algorithm is also easy to understand.
"""
function kinship(ped, i::Int, j::Int, dic::Dict{iPair, Float64})
    (i == 0 || j == 0) && return 0
    ip, im = ped[i, :]
    if i == j
        haskey(dic, (i, i)) || (dic[(i, i)] = 1 + .5kinship(ped, ip, im, dic))
        return dic[(i, i)]
    end
    if i < j
        jp, jm = ped[j, :]
        haskey(dic, (i, jp)) || (dic[(i, jp)] = kinship(ped, i, jp, dic))
        haskey(dic, (i, jm)) || (dic[(i, jm)] = kinship(ped, i, jm, dic))
        return .5(dic[(i, jp)] + dic[(i, jm)])
    end
    return .5(kinship(ped, j, ip, dic) + kinship(ped, j, im, dic))
end

#=
"""
    function ped_F(ped; force = false, mdc = 1_000_000)
- When column `:F` is not in DataFrame `ped`, this function calculate
inbreeding coefficients of every ID, and add a `:F` column in `ped`.
- If `:F` exists, and `force = true`, drop `:F` from `ped` and do above.
- Or do nothing.

The relationship dictionary is limited to 1M to save memory
when the pedigree is too large.
"""
function ped_F(ped; force = false, mdc = 1_000_000)
    if "F" ∈ names(ped)
        force ? select!(ped, Not([:F])) : return
    end
    N = nrow(ped)
    F = zeros(N)
    dic = Dict{Tuple{Int, Int}, Float64}()
    @showprogress for i in 1:N
        F[i] = kinship(ped, i, i, dic) - 1
    end
    ped.F = F
    nothing
end

"""
    function ped_D(ped; force = false)
Calculate diagonals of `D` for `A` or `A⁻¹` calculation.
It is computed as a vector and a new column of `ped` of name `:D`.
"""
function ped_D(ped; force = false)
    if "D" ∈ names(ped)
        force ? select!(ped, Not([:D])) : return
    end
    "F" ∈ names(ped) || ped_F(ped)
    N = nrow(ped)
    D = .5ones(N)
    @showprogress for i in 1:N
        pa, ma = ped[i, :]
        vp = (pa == 0) ? -.25 : .25ped.F[pa]
        vm = (ma == 0) ? -.25 : .25ped.F[ma]
        D[i] -= (vp + vm)
    end
    ped.D = D
    nothing
end

"""
    function _D4A_ni(ped; inverse = false)
Construct a `D` diagonal matrix for `A` related calculation.
This one is for test only, and is hence internal.
"""
function _D4A_ni(ped; inverse = false)
    D = ones(nrow(ped))
    for (i, (pa, ma, ..)) in enumerate(eachrow(ped))
        pa > 0 && (D[i] -= .25)
        ma > 0 && (D[i] -= .25)
    end
    inverse && (D = 1 ./ D)
    Diagnal(D)
end

"""
    function _pushRCV!(R, C, V, r, c, v)
Push a value row, column, and value into vector `R`, `C` and `V`.
The vectors are later to be used to construct a sparse matrix.
"""
function _pushRCV!(R, C, V, r, c, v)
    push!(R, r)
    push!(C, c)
    push!(V, v)
end

"""
    function T4AI(ped)
Give a pedigree DataFrame, with its first 2 column as `pa`, and `ma`,
this function return the T matrix used for A⁻¹ calculation.
"""
function T4AI(ped)
    N = nrow(ped)
    # R, C, V: row, column and value specifying values in a sparse matrix
    # 3N are enough, as diagonal → N, all parents known → another 2N.
    R, C, V = Int[], Int[], Float64[]
    for (id, (pa, ma, ..)) in enumerate(eachrow(ped))
        pa > 0 && _pushRCV!(R, C, V, id, pa, -.5)
        ma > 0 && _pushRCV!(R, C, V, id, ma, -.5)
        _pushRCV!(R, C, V, id, id, 1.)
    end
    T = sparse(R, C, V)
end

"""
    function T4A(ped; m = 1000)
Calculate `T`, which can be used to construct `A` as `TDT'`.
`T` can also be deemed as genetic contribution matrix.
"""
function T4A(ped; m = 1000)
    N = nrow(ped)
    N > m && error("Pedigree size ($N > $m), too big")
    T = zeros(N, N) + I(N)
    for (i, (pa, ma, ..)) in enumerate(eachrow(ped))
        if pa > 0
            for j in 1:i-1
                T[i, j] += .5T[pa, j]
            end
        end
        if ma > 0
            for j in 1:i-1
                T[i, j] += .5T[ma, j]
            end
        end
    end
    T     
end

"""
    function Amat(ped; m = 1000)
Given a pedigree `ped`,
this function returns a full numerical relationship matrix, `A`.
This function is better for small pedigrees, and for demonstration only.
The maximal matrix size is thus limited to 1000.
One can try to set `m` to a bigger value if RAM is enough.
"""
function Amat(ped; m = 1000)
    N = nrow(ped)
    N > m && error("Pedigree size ($N > $m) too big")
    A = zeros(N, N) + I(N)
    for (i, (pa, ma, ..)) in enumerate(eachrow(ped))
        pa * ma ≠ 0 && (A[i, i] += .5A[pa, ma])
        for j in 1:i-1
            pa ≠ 0 && (A[i, j]  = 0.5A[j, pa])
            ma ≠ 0 && (A[i, j] += 0.5A[j, ma])
            A[j, i] = A[i, j]
        end
    end
    A
end
        
"""
    function A22(ped, idc, list::Vector{String}, oo::String)
Construct sub matrix of `A` of ID in `list`.
Dictionary `idc` records the coded ID, i.e., row numbe, in `ped`.
It was created when reading `ped`.
To prevent out of memory error, results are written to file `oo`.
"""
function A22(ped, idc, list::Vector{String}, oo::String)
    N = length(list)
    ilst = zeros(Int, N)
    i = 0
    for id in list
        i += 1
        haskey(idc, id) ? (ilst[i] = idc[id]) : error("$id not found in idc")
    end
    df = DataFrame(x = Int[], y = Int[], v = Float64[])
    th = Threads.nthreads() * 10000 # number of pairs to calculate per batch
    open(oo, "w+") do io
        write(io, [N, N, Fio.typec(Float64)])
        rst = Mmap.mmap(io, Matrix{Float64}, (N, N), 24)
        @showprogress 1 "Parallel computing A[i, j]" for i in 1:N
            for j in 1:i
                push!(df, (ilst[i], ilst[j], 0))
                if nrow(df) == th
                    Threads.@threads for k in 1:th
                        df.v[k] = kinship(ped, df.x[k], df.y[k])
                    end
                    for (x, y, v) in eachrow(df)
                        A[x, y] = v
                    end
                    empty!(df)
                end
            end
        end
        @info "Finalizing ..."
        if !isempty(pairs)
            Threads.@threads for k in 1:nrow(df)
                df.v[k] = kinship(ped, df.x[k], df.y[k])
            end
            for (x, y, v) in eachrow(df)
                A[x, y] = v
            end
        end
    end
    nothing
end

"""
    function Ai(ped)
Given a pedigree `ped`, this function return a sparse matrix of `A⁻¹`,
where `A` is the numerical relationship matrix.
"""
function Ai(ped)
    T = T4AI(ped)
    "D" ∈ names(ped) || ped_D(ped)
    D = Diagonal(1. ./ ped.D)
    T'D*T
end

"""

"""
function effID(ped, idc, lst)
    epp = Set{Int}()
    id = Set{Int}()
    for s in lst
        push!(id, idc[s])
    end
    ng = 1
    while !isempty(id)
        union!(epp, unique(id))
        pm = Set{Int}()
        ng += 1
        for i in id
            ped.pa[i] == 0 || push!(pm, ped.pa[i])
            ped.ma[i] == 0 || push!(pm, ped.ma[i])
        end
        id = pm
        @show ng, length(epp), length(pm)
    end
    sort(collect(epp))
end
=#