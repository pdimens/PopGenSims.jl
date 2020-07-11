export cross, backcross


function haploid_cross!(data::DataFrame, parent::String, parent2::String; n::Int = 100)
    iter_df = DataFrames.groupby(data, :name)
    for simulation in iter_df
        all_alleles = getindex.(collect.(zip.(p1,p2)), 1)
        simulation.genotype .= Tuple(
                                Base.Iterators.flatten(rand(all_alleles, 1)) |> collect
                            )
    end
    return data
end

function polyploid_cross!(data::DataFrame, p1::T, p2::T; n::Int = 100, ploidy::Int = 2) where T <: GenoArray
    iter_df = DataFrames.groupby(data, :name)
    n_alleles = ploidy ÷ 2
    for simulation in iter_df
        simulation.genotype .= Tuple.(
                                zip(
                                    Base.Iterators.flatten(rand.(p1, n_alleles)) |> collect,
                                    Base.Iterators.flatten(rand.(p2, n_alleles)) |> collect
                                )
                            )
    end
    return data
end

"""
    cross(data::PopData, parent::String, parent2::String; n::Int = 100, ploidy::Union{Int, String}= 2)

Simulate a breeding cross between individuals `parent` and `parent2` from a `PopData` object.
Returns PopData consisting of `n` offspring resulting from the cross.
#### Keyword Arguments
- `n` : Integer of number of offspring to generate (default: `100`)
- `ploidy`: Integer or String of the ploidy of the samples. 
"""
function cross(data::PopData, parent::String, parent2::String; n::Int = 100)
     # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)
    # check for parents not having mixed ploidy
    length(unique(length.(skipmissing(p1)))) != 1 && error("Parent $parent has mixed ploidy, which is unsupported")
    length(unique(length.(skipmissing(p2)))) != 1 && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # Get the ploidy value & check for equal ploidy
    p1_ploidy = length.(skipmissing(p1))[1]
    p2_ploidy = length.(skipmissing(p2))[1]
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    loci = levels(data.loci.locus)
    
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    parents = Vector{Tuple{String, String}}(undef, length(out_loci))
    parents .= Ref((parent, parent2))
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :parents => parents, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :parents, :locus], compress = true)

    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent, parent2, n = n)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_df, p1, p2, n = n, ploidy = p1_ploidy)
    else
        error("Currently supported ploidy: 1, 2, 4, 6, 8")
    end
end


"""
    backcross(parent_1_data::PopData, parent_2_data::DataFrame; parent1::String, parent2::String, n::Int = 100)
Cross a parent with an offspring generated from using `cross`.
"""
function backcross(parent_1_data::PopData, parent_2_data::DataFrame; parent1::String, parent2::String, n::Int = 100)
    p1 = get_genotypes(parent_1_data, parent1)
    p2 = @view parent_2_data[parent_2_data.name .== parent2, :genotype]

    # check for parents not having mixed ploidy
    length(unique(length.(skipmissing(p1)))) != 1 && error("Parent $parent has mixed ploidy, which is unsupported")
    length(unique(length.(skipmissing(p2)))) != 1 && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # Get the ploidy value & check for equal ploidy
    p1_ploidy = length.(skipmissing(p1))[1]
    p2_ploidy = length.(skipmissing(p2))[1]
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    loci = levels(parent_1_data.loci.locus)
    
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    parents = Vector{Tuple{String, String}}(undef, length(out_loci))
    parents .= Ref((parent1, parent2))
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :parents => parents, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :parents, :locus], compress = true)

    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent, parent2, n = n)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_df, p1, p2, n = n, ploidy = p1_ploidy)
    else
        error("Currently supported ploidy: 1, 2, 4, 6, 8")
    end
end