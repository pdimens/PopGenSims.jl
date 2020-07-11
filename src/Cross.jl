export cross


"""

"""
function haploid_cross(data::PopData, parent::String, parent2::String; n::Int = 100)
    # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)
    loci = levels(data.loci.locus)

    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :locus], compress = true)
    for simulation in iter_df
        all_alleles = getindex.(collect.(zip.(p1,p2)), 1)
        simulation.genotype .= Tuple(
                                Base.Iterators.flatten(rand(all_alleles, 1)) |> collect
                            )
    end
end

function diploid_cross(data::PopData, parent::String, parent2::String; n::Int = 100)
    # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)

    loci = levels(data.loci.locus)
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :locus], compress = true)
    
    iter_df = DataFrames.groupby(out_df, :name)
    for simulation in iter_df
        simulation.genotype .= Tuple.(
                                zip(
                                    Base.Iterators.flatten(rand.(p1, 1)) |> collect,
                                    Base.Iterators.flatten(rand.(p2, 1)) |> collect
                                )
                            )
    end
    return out_df
end

function tetraploid_cross(data::PopData, parent::String, parent2::String; n::Int = 100)
    # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)

    loci = levels(data.loci.locus)
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :locus], compress = true)
    
    iter_df = DataFrames.groupby(out_df, :name)
    for simulation in iter_df
        simulation.genotype .= Tuple.(
                                zip(
                                    Base.Iterators.flatten(rand.(p1, 2)) |> collect,
                                    Base.Iterators.flatten(rand.(p2, 2)) |> collect
                                )
                            )
    end
    return out_df
end

function hexaploid_cross(data::PopData, parent::String, parent2::String; n::Int = 100)
    # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)

    loci = levels(data.loci.locus)
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :locus], compress = true)
    
    iter_df = DataFrames.groupby(out_df, :name)
    for simulation in iter_df
        simulation.genotype .= Tuple.(
                                zip(
                                    Base.Iterators.flatten(rand.(p1, 3)) |> collect,
                                    Base.Iterators.flatten(rand.(p2, 3)) |> collect
                                )
                            )
    end
    return out_df
end

function octaploid_cross(data::PopData, parent::String, parent2::String; n::Int = 100)
    # get parental genotypes
    p1 = get_genotypes(data, parent)
    p2 = get_genotypes(data, parent2)

    loci = levels(data.loci.locus)
    # pre-allocate all output information
    out_loci = fill(loci, n) |> Base.Iterators.flatten |> collect
    out_offspring = fill.(["offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_geno = similar(p1, length(loci) * n)
    out_df = DataFrame(:name => out_offspring, :locus => out_loci, :genotype => out_geno)
    categorical!(out_df, [:name, :locus], compress = true)
    
    iter_df = DataFrames.groupby(out_df, :name)
    for simulation in iter_df
        simulation.genotype .= Tuple.(
                                zip(
                                    Base.Iterators.flatten(rand.(p1, 4)) |> collect,
                                    Base.Iterators.flatten(rand.(p2, 4)) |> collect
                                )
                            )
    end
    return out_df
end

function cross(data::PopData, parent::String, parent2::String; n::Int = 100, ploidy::Union{Int, String}= 2)

end