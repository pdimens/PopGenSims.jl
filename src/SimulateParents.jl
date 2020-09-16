#TODO in PopGen.jl, rename alleles => allele_pool and replace with this faster method
#=
function allele_pool(locus::T) where T <: GenoArray
    reduce(vcat, collect.(skipmissing(locus)))
end
=#
function allele_pool(locus::T) where T <: GenoArray
    Tuple(Base.Iterators.flatten(skipmissing(locus)))
end

function allele_pool(data::PopData)
    # index dataframe by locus
    idx_df = groupby(data.loci, [:locus])
    # get the data type
    markertype = skipmissing(idx_df[1].genotype) |> first |> eltype
    # instantiate dict to store alleles
    allele_dict = Dict{String,NTuple}()
    # pull out loci names
    loc = getindex.(keys(idx_df), :locus)
    [allele_dict[i] = allele_pool(idx_df[(;locus = i)].genotype) for i in loc]
    return String.(loc), allele_dict
end

function simulate_parent(alleles::Dict{String,NTuple}, loc::Vector{String}; ploidy::Int = 2)
    map(i -> rand(alleles[i], ploidy) ,loc)
end

function cross(parent1::Vector{Vector{T}}, parent2::Vector{Vector{T}}) where T <: Signed
    p1_contrib = rand.(parent1)
    p2_contrib = rand.(parent2)
    sort.(zip(p1_contrib, p2_contrib))
end

#= wont need
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end
=#
