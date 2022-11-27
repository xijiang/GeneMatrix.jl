# © 2022, Xijiang Yu.  MIT license
# Numerical Relationship Matrix

"""
    nrm(ped)

Calculate the numerical relationship matrix of the given pedigree `ped`.
The pedigree should be a `DataFrame` with its first 2 columns as `Pa` 
and `Ma` of integers. The row numbers are deemed as `ID` numbers.  This
pedigree must be sorted, such that an `ID`, or row numbers always appear
after `Pa` and `Ma`.

An sorted pedigree that meet above can be obtained with function `sortped`.
"""
function nrm(ped)
end

function nrminv(ped)
    
end

function sortped(file)
end
