module PopGenSims

using DataFrames, PooledArrays, StaticArrays
using StatsBase: sample, Weights
using PopGen:
    allele_freq,
    Genotype,
    GenoArray,
    get_genotypes,
    PopData,
    read_from,
    sort,
    write_to

include("Cross.jl")
include("Samples.jl")
include("Sibship.jl")
include("Utils.jl")

end
