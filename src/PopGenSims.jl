module PopGenSims

using DataFrames, PooledArrays, StaticArrays
using StatsBase: sample, Weights
using RandomNumbers.Xorshifts
using PopGen:
    PopData,
    Genotype,
    GenoArray,
    get_genotypes,
    read_from,
    write_to,
    sort

include("Cross.jl")
include("Samples.jl")
include("Sibship.jl")
include("Utils.jl")

end
