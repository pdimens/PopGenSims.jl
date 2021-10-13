module PopGenSims

using DataFrames, PooledArrays, StaticArrays
using StatsBase: sample, Weights
import PopGenCore: read, write
using PopGenCore:
    Genotype,
    GenoArray,
    get_genotypes,
    PopData, 
    sort

using PopGen: allele_freq

include("Cross.jl")
include("Samples.jl")
include("Sibship.jl")
include("Utils.jl")

end
