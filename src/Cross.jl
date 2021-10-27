export cross, backcross

function sample_genotype(geno::T, n_alleles::Int) where T<:Genotype
    sample(SVector(geno), n_alleles, replace = false)
end

function sample_genotype(geno::Missing, n_alleles::Int)
    return [missing]
end

function haploid_cross!(data::DataFrame, p1::T, p2::T) where T <: GenoArray
    iter_df = DataFrames.groupby(data, :name)
    for simulation in iter_df
        all_alleles = getindex.(collect.(zip.(p1,p2)), 1)
        offspring_geno = Tuple(Base.Iterators.flatten(rand.(all_alleles, 1)) |> collect)
        miss_idx = reduce(union, findall.(i -> i == (missing,), offspring_geno))
        simulation.genotype[Not(miss_idx)] .= offspring_geno[Not(miss_idx)]
    end
    return data
end

function polyploid_cross!(data::DataFrame, p1::T, p2::T; ploidy::Signed) where T <: GenoArray
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
    err = ""
    err *= parent1 ∉ (@view data.sampleinfo[!, :name]) ? "$parent1 " : ""
    err *= parent2 ∉ (@view data.sampleinfo[!, :name]) ? "$parent2" : ""
    err != "" && error("One or more parents not found in PopData: $err")
    
    # Get the ploidy value & check for equal ploidy
    p1_ploidy = data.sampleinfo.ploidy[data.sampleinfo.name .== parent1] |> first
    p2_ploidy = data.sampleinfo.ploidy[data.sampleinfo.name .== parent2] |> first
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    # check for parents not having mixed ploidy
    p1_ploidy isa AbstractVector && error("Parent $parent1 has mixed ploidy, which is unsupported")
    p2_ploidy isa AbstractVector && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # get parental genotypes
    p1 = genotypes(data, parent1)
    p2 = genotypes(data, parent2)

    loci = data.locusinfo.locus
    
    # pre-allocate all output information
    out_loci_names = repeat(loci, outer = n)
    _padding = length(string(n))
    out_offspring = repeat(["$generation" * "_" * lpad("$i", _padding, "0") for i in 1:n], inner = length(loci))
    out_population = fill(generation, n * length(loci))
    out_geno = similar(p1, n * length(loci))
    out_loci = DataFrame(
        :name => PooledArray(out_offspring, compress = true), 
        :population => PooledArray(out_population, compress = true), 
        :locus => PooledArray(out_loci_names, compress = true), 
        :genotype => out_geno
        )
    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent1, parent2)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_loci, p1, p2, ploidy = p1_ploidy)
    else
        throw(MethodError("Currently supported ploidy: 1, 2, 4, 6, 8"))
    end
    out = PopData(out_loci)
    insertcols!(out.sampleinfo, :parents => PooledArray(fill((parent1,parent2), n), compress = true))
    return out
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
- `generation` : A string to assign `population` identity and name prefix to the offspring (default: `"F1"`)
"""
function cross(parent_1::Pair, parent_2::Pair; n::Int = 100, generation::String = "F1")
    parent_1_data = parent_1.first
    parent_2_data = parent_2.first

    loci = parent_1_data.locusinfo.locus
    loci_p2 = parent_2_data.locusinfo.locus
    length(loci) != length(loci_p2) && error("Both datasets must have the same number of loci. $parent1 : $length(loci) | $parent2 : $length(loci_p2")
    
    # verify identical loci
    loci_check = loci .!= loci_p2
    culprits_p1 = loci[loci_check]
    culprits_p2 = loci_p2[loci_check]
    culp_print = "Parent 1\tParent 2" * "\n---------\t---------\n" * join("$i\t$j\n" for (i,j) in zip(culprits_p1, culprits_p2))
    length(culprits_p1) > 0 && error("Both datasets must have loci in the same order. Loci causing this error:\n" * culp_print)

    parent1 = parent_1.second
    parent2 = parent_2.second

    # check for presence of parents
    parent1 ∉ (@view parent_1_data.sampleinfo[!, :name]) && error("$parent1 not found in PopData")
    parent2 ∉ (@view parent_2_data.sampleinfo[!, :name]) && error("$parent2 not found in PopData")

    # Get the ploidy value & check for equal ploidy
    p1_ploidy = parent_1_data.sampleinfo.ploidy[parent_1_data.sampleinfo.name .== parent1] |> first
    p2_ploidy = parent_2_data.sampleinfo.ploidy[parent_2_data.sampleinfo.name .== parent2] |> first
    p1_ploidy != p2_ploidy && error("Parents must have identical ploidy. Parent1 = $p1_ploidy | Parent2 = $p2_ploidy")

    # check for parents not having mixed ploidy
    p1_ploidy isa AbstractVector && error("Parent $parent1 has mixed ploidy, which is unsupported")
    p2_ploidy isa AbstractVector && error("Parent $parent2 has mixed ploidy, which is unsupported")

    # get parental genotypes
    p1 = genotypes(parent_1_data, parent1)
    p2 = genotypes(parent_2_data, parent2)

    # pre-allocate all output information
    out_loci_names = repeat(loci, outer = n)
    _padding = length(string(n))
    out_offspring = repeat(["$generation" * "_" * lpad("$i", _padding, "0") for i in 1:n], inner = length(loci))
    out_population = fill(generation, n * length(loci))
    out_geno = similar(p1, n * length(loci))
    out_loci = DataFrame(:name => out_offspring, :population => out_population, :locus => out_loci_names, :genotype => out_geno)
    out_loci.name = PooledArray(out_loci.name, compress = true)
    out_loci.population = PooledArray(out_loci.population, compress = true)
    out_loci.locus = PooledArray(out_loci.locus, compress = true)

    # perform the cross
    if p1_ploidy == 1 
        haploid_cross!(data, parent1, parent2)
    elseif p1_ploidy ∈  [2, 4, 6, 8] 
        polyploid_cross!(out_loci, p1, p2, ploidy = p1_ploidy)
    else
        error("Currently supported ploidy: 1, 2, 4, 6, 8")
    end
    out = PopData(out_loci)
    insertcols!(out.sampleinfo, :parents => fill((parent1,parent2), n))
    return out
end