module GeneticCrosses

using DataFrames, PooledArrays
using PopGen:
    PopData,
    Genotype,
    GenoArray,
    get_genotypes,
    read_from,
    write_to,
    nancycats

include("Cross.jl")

end

