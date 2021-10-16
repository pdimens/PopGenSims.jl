module PopGenSims

using DataFrames, PooledArrays, StaticArrays
using StatsBase: sample, Weights
import PopGenCore: read, write
using PopGenCore:
    allele_freq,
    copy,
    Genotype,
    GenoArray,
    get_genotypes,
    PopData, 
    sort

include("Cross.jl")
include("Samples.jl")
include("Sibship.jl")
include("Utils.jl")

end
