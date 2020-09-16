module GeneticCrosses

using DataFrames, PooledArrays, StaticArrays
using PopGen:
    PopData,
    Genotype,
    GenoArray,
    get_genotypes,
    read_from,
    write_to,
    nancycats, 
    sort

include("Cross.jl")

end

