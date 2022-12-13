# the two functions below return a dictionary with 81 items.
# they are of opposite directions.

"""
    bit2floatDict()

This function uses the plink coding convension, except it doesn't accept
missing values.
The returned dictionary translate a byte into 4 Float64.
"""
function bit2floatDict()
    dic = Dict{Int8, NTuple{4, Float64}}()
    v4 = [0, 0, 0, 0]
    wd = Int8.([0, 2, 3]) # 00->0, 10->1, 11->2
    for _ in 1:81
        byte = Int8(0)
        ns = 0  # number of shifts
        for i in v4
            gt = wd[i+1]
            gt <<= ns
            byte |= gt
            ns += 2
        end
        dic[byte] = Tuple(Float64.(v4))
        v4[1] += 1
        for i in 1:3
            v4[i+1] += v4[i] ÷ 3
            v4[i] %= 3
        end
    end
    dic
end

"""
    char2bitDict()

This function uses the plink coding convension, except it doesn't accept
missing values.
The returned dictionary codes 4 Int8 genotypes into a byte.
"""
function char2bitDict()
    dic = Dict{NTuple{4, Int8}, Int8}()
    v4 = [0, 0, 0, 0]
    wd = Int8.([0, 2, 3]) # 00->0, 10->1, 11->2
    for _ in 1:81
        byte = Int8(0)
        ns = 0  # number of shifts
        for i in v4
            gt = wd[i+1]
            gt <<= ns
            byte |= gt
            ns += 2
        end
        dic[Tuple(Float64.(v4))] = byte
        v4[1] += 1
        for i in 1:3
            v4[i+1] += v4[i] ÷ 3
            v4[i] %= 3
        end
    end
    dic
end

"""
    char2bit(mat::Matrix{Int8})

Code genotype matrix `mat` int plink bed type matrix.
This function always translate the genotypes by columns.
Missing genotype is not allowed.
"""
function char2bit(mat::Matrix{Int8})
    r, c = size(mat)
    m, n = divrem(r, 4)
    bat = n>0 ? zeros(Int8, m+1, c) : zeros(Int8, m, c)
    dic = char2bitDict()
    for j in 1:c
        for i in 1:m
            t4 = Tuple(mat[4i-3:4i, j])
            bat[i, j] = dic[t4]
        end
    end
    if n > 0
        a = zeros(4 - n)
        l = 4m+1:4m+n
        for j in 1:c
            t4 = Tuple(vcat(mat[l, j], a))
            bat[end, j] = dic[t4]
        end
    end
    bat
end

"""
    bit2float(mat::Matrix{Int8}, nrow)
Convert plink 2bit file into Float64.
"""
function bit2float(mat::Matrix{Int8}, nrow)
    r, c = size(mat)
    4r-4 < nrow ≤ 4r || error("Row number $nrow error.")
    fat = zeros(nrow, c)
    dic = bit2floatDict()
    m = 4r > nrow ? r - 1 : r
    for i in 1:m
        for j in 1:c
            copyto!(view(fat, 4i-3:4i, j), dic[mat[i, j]])
        end
    end
    if 4r > nrow
        x = 1:(nrow+4-4r)
        y = (4r - 3):nrow
        for j in 1:c
            copyto!(view(fat, y, j), dic[mat[end, j]][x])
        end
    end
    fat
end
