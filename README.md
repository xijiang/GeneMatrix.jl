# GeneMatrix

This package is iniatially contended for the constructions various
matrices used in animal breeding and genetics.  It also provides a set
of file I/O functions to store and retrieve the matrices to and from a
storage, e.g, HDD.  The reason for the latter is that saving the
matrix might occupy too much disk space, e.g., the 3-column,
row-column-value storage.  It also prevents precision loss while store
float number in text format.

## File I/O

The storage format is called `XY` format.  Each contains a header, as
defined by ``struct xyheader``.  The starting 3 characters of a header
are "XY ".  The 4th character is any one of "FLUS", that is `F`ull
matrix, `L`ower triangle, `U`pper triangle, or lower triangle of a
`S`ymmetric matrix.  the 5th character specify if the genotypes matrix
is of `nLoci $\times$ nID` (`N`), or `T`ransposed if applies.  The 6th
character specify element types of the matrix.  The currently
supported types are ``Bool, Int8, Int16, Int32, Int64, Int128, UInt8,
UInt16, UInt32, UInt64, UInt128, Float16, Float32, Float64``, though
some of them may never be used.  If a new type is to be added, it must
appended to the end of above.

### I/O functions

The function for reading is `readxy`.  It only take one argument of
file name.  The one for writing is `writexy(file, mat; mattp = 'F',
trans = 'N')`.  It write matrix `mat` into `file`.  If `mat` is
symmetric, then `mattp` is overruled as `S`.

The package also provides two other I/O functions to input and output
plink bed files.

- *To be added*

## Matrix construction

### Numerical relationship matrix related

### Genome relationship matrix related

## Pedigree

## Compact genotype storage and manipulation methods

## Todo

- `nrm`: numerical relationship matrix
- `nrminv`: the inverse of nrm
- `sortped`, recode pedigree to 1:nID.
- `readmacs`, read simulation results from MaCS.
- `bed2int8`, convert bed file formats to Int8
- tests merging.
