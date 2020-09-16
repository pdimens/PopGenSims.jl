export cross, backcross

function sample_genotype(geno::T, n_alleles::Int) where T<:Genotype
    sample(geno, n_alleles, replace = false)
end

function sample_genotype(geno::Missing, n_alleles::Int)
    return [missing]
end

function haploid_cross!(data::DataFrame, p1::T, p2::T; n::Int) where T <: GenoArray
    iter_df = DataFrames.groupby(data, :name)
    for simulation in iter_df
        all_alleles = getindex.(collect.(zip.(p1,p2)), 1)
        offspring_geno = Tuple(Base.Iterators.flatten(rand.(all_alleles, 1)) |> collect)
        miss_idx = reduce(union, findall.(i -> i == (missing,), offspring_geno))
        simulation.genotype[Not(miss_idx)] .= offspring_geno[Not(miss_idx)]
    end
    return data
end

function polyploid_cross!(data::DataFrame, p1::T, p2::T; n::Int, ploidy::Int) where T <: GenoArray
    iter_df = DataFrames.groupby(data, :name)
    n_alleles = ploidy ÷ 2
    for simulation in iter_df
        p1_contrib = Base.Iterators.flatten(sample_genotype.(p1, n_alleles)) |> collect
        p2_contrib = Base.Iterators.flatten(sample_genotype.(p2, n_alleles)) |> collect
        miss_idx = reduce(union, findall.(ismissing, [p1_contrib, p2_contrib]))
        offspring_geno = zip(p1_contrib, p2_contrib) |> collect
        simulation.genotype[Not(miss_idx)] .= offspring_geno[Not(miss_idx)]
    end
    return data
end

"""
    cross(data::PopData, parent1::String, parent2::String; n::Int = 100, generation::String = "F1")

Simulate a breeding cross between individuals `parent1` and `parent2` from a `PopData` object.
Returns PopData consisting of `n` offspring resulting from the cross.
#### Keyword Arguments
- `n` : Integer of number of offspring to generate (default: `100`)
- `generation` : A string to assign `population` identity to the offspring (default: `"F1"`)
"""
function cross(data::PopData, parent1::String, parent2::String; n::Int = 100, generation::String = "F1")
    # check for presence of parents
    parent1 ∉ (@view data.meta[!, :name]) && error("$parent not found in PopData")
    parent2 ∉ (@view data.meta[!, :name]) && error("$parent2 not found in PopData")
    
    # get parental genotypes
    p1 = get_genotypes(data, parent1)
    p2 = get_genotypes(data, parent2)

    # check for parents not having mixed ploidy
    length(unique(length.(skipmissing(p1)))) != 1 && error("Parent $parent has mixed ploidy, which is unsupported")
    length(unique(length.(skipmissing(p2)))) != 1 && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # Get the ploidy value & check for equal ploidy
    p1_ploidy = length.(skipmissing(p1)) |> first
    p2_ploidy = length.(skipmissing(p2)) |> first
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    loci = unique(data.loci.locus)
    
    # pre-allocate all output information
    out_loci_names = fill(loci, n) |> Base.Iterators.flatten |> collect
    #parents = Vector{Tuple{String, String}}(undef, n)

    out_offspring = fill.(["$generation" * "_offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_population = fill(generation, n * length(loci))
    out_geno = similar(p1, n * length(loci))
    out_loci = DataFrame(:name => out_offspring, :population => out_population, :locus => out_loci_names, :genotype => out_geno)
    #categorical!(out_loci, [:name, :population, :locus], compress = true)
    out_loci.name = PooledArray(out_loci.name)
    out_loci.population = PooledArray(out_loci.population)
    out_loci.locus = PooledArray(out_loci.locus)
    
    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent, parent2, n = n)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_loci, p1, p2, n = n, ploidy = p1_ploidy)
    else
        error("Currently supported ploidy: 1, 2, 4, 6, 8")
    end
    out_meta = DataFrame(
        :name => unique(out_loci.name),
        :ploidy => fill(p1_ploidy, n),
        :population => fill(generation, n),
        :latitude => Vector{Union{Missing, Float32}}(undef, n),
        :longitude => Vector{Union{Missing, Float32}}(undef, n),
        :parents => fill((parent,parent2), n)
    )
    PopData(out_meta, out_loci)
end


"""
    cross(parent_1::Pair, parent_2::Pair, n::Int = 100, generation::String = "F1")

Simulate a breeding cross between individuals `parent` and `parent2` from two different `PopData` objects.
Returns PopData consisting of `n` offspring resulting from the cross. `parent_1_data` and `parent_2_data` 
are positional arguments, therefore they must be written without keywords and in the order of parents 1, parent 2. 
#### Keyword Arguments
- `parent_1` : Pair of `PopData => "Parent1Name"`
- `parent_2` : Pair of `PopData => "Parent1Name"`
- `n` : Integer of number of offspring to generate (default: `100`)
- `generation` : A string to assign `population` identity to the offspring (default: `"F1"`)
"""
function cross(parent_1::Pair, parent_2::Pair; n::Int = 100, generation::String = "F1")
    parent_1_data = parent_1.first
    parent_2_data = parent_2.first

    parent1 = parent_1.second
    parent2 = parent_2.second

    # check for presence of parents
    parent1 ∉ (@view parent_1_data.meta[!, :name]) && error("$parent not found in PopData")
    parent2 ∉ (@view parent_2_data.meta[!, :name]) && error("$parent2 not found in PopData")
    
    # get parental genotypes
    p1 = get_genotypes(parent_1_data, parent1)
    p2 = get_genotypes(parent_2_data, parent2)

    # check for parents not having mixed ploidy
    length(unique(length.(skipmissing(p1)))) != 1 && error("Parent $parent has mixed ploidy, which is unsupported")
    length(unique(length.(skipmissing(p2)))) != 1 && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # Get the ploidy value & check for equal ploidy
    p1_ploidy = length.(skipmissing(p1)) |> first
    p2_ploidy = length.(skipmissing(p2)) |> first
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    # verify identical loci
    loci = unique(parent_1_data.loci.locus)
    loci_p2 = unique(parent_2_data.loci.locus)
    length(loci) != length(loci_p2) && error("Both parents must have the same number of loci. $parent1 : $length(loci) | $parent2 : $length(loci_p2")
    loci_check = loci .!= loci_p2
    culprits_p1 = loci[loci_check]
    culprits_p2 = loci_p2[loci_check]
    culp_print = "Parent 1\tParent 2" * "\n---------\t---------\n" * join("$i\t$j\n" for (i,j) in zip(culprits_p1, culprits_p2))
    length(culprits_p1) > 0 && error("Both datasets must have loci in the same order. Loci causing this error:\n" * culp_print)
    
    # pre-allocate all output information
    out_loci_names = fill(loci, n) |> Base.Iterators.flatten |> collect
    #parents = Vector{Tuple{String, String}}(undef, n)

    out_offspring = fill.(["$generation" * "_offspring_$i" for i in 001:n], length(loci)) |> Base.Iterators.flatten |> collect
    out_population = fill(generation, n * length(loci))
    out_geno = similar(p1, n * length(loci))
    out_loci = DataFrame(:name => out_offspring, :population => out_population, :locus => out_loci_names, :genotype => out_geno)
    #categorical!(out_loci, [:name, :population, :locus], compress = true)
    out_loci.name = PooledArray(out_loci.name)
    out_loci.population = PooledArray(out_loci.population)
    out_loci.locus = PooledArray(out_loci.locus)

    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent, parent2, n = n)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_loci, p1, p2, n = n, ploidy = p1_ploidy)
    else
        error("Currently supported ploidy: 1, 2, 4, 6, 8")
    end
    out_meta = DataFrame(
        :name => unique(out_loci.name),
        :ploidy => fill(p1_ploidy, n),
        :population => fill(generation, n),
        :latitude => Vector{Union{Missing, Float32}}(undef, n),
        :longitude => Vector{Union{Missing, Float32}}(undef, n),
        :parents => fill((parent1,parent2), n)
    )
    PopData(out_meta, out_loci)
end