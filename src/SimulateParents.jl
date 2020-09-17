#TODO in PopGen.jl, rename alleles => allele_pool and replace with this faster method
#=
function allele_pool(locus::T) where T <: GenoArray
    reduce(vcat, collect.(skipmissing(locus)))
end
=#
function allele_pool(locus::T) where T <: GenoArray
    Tuple(Base.Iterators.flatten(skipmissing(locus)))
end

#= wont need
function Base.sort(x::NTuple{N,T}) where N where T <: Signed 
    Tuple(sort(SVector(x)))
end
=#

function allele_pool(data::PopData)
    # index dataframe by locus
    idx_df = groupby(data.loci, [:locus])
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
    ploidy = length(first(parent1))
    ploidy == 1 && error("Haploid crosses are not yet supported. Please file and issue or pull request")
    if ploidy == 2
        p1_contrib = rand.(parent1)
        p2_contrib = rand.(parent2)
        geno_out = sort.(zip(p1_contrib, p2_contrib))
    elseif iseven(ploidy)
        n_allele = ploidy ÷ 2
        p1_contrib = sample.(parent1, n_allele, replace = false)
        p2_contrib = sample.(parent2, n_allele, replace = false)
        geno_out = Tuple.(sort!.(append!.(p1_contrib, p2_contrib)))
    else
        # special method to provide a 50% chance of one parent giving more alleles than the other
        rng = rand()
        contrib_1 = ploidy ÷ 2
        contrib_2 = ploidy - contrib_1
        p1_contrib = rng > 0.5 ? sample.(parent1, contrib_1, replace = false) : sample.(parent1, contrib_2, replace = false)
        p2_contrib = rng > 0.5 ? sample.(parent2, contrib_2, replace = false) : sample.(parent2, contrib_1, replace = false)
        geno_out = Tuple.(sort!.(append!.(p1_contrib, p2_contrib)))
    end
    return geno_out
end

function parentoffspring(data::PopData; n::Int = 100, ploidy::Int = 2)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_parent(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_parent") => cross(p1, p2))
        insertcols!(out_df, Symbol(prefix * "_offspring") => Tuple.(p1))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "parent_offspring")
    out_df.name = PooledArray(out_df.name)
    out_df.population = PooledArray(out_df.population)
    out_df.locus = PooledArray(out_df.locus)
    
    meta_df = DataFrame(
        :name => unique(out_df.name), 
        :population => "parent_offspring", 
        :ploidy => ploidy, 
        :longitude => Vector{Union{Missing, Float32}}(undef, n*2),
        :latitude => Vector{Union{Missing, Float32}}(undef, n*2),
        )
    return PopData(meta_df, out_df)
end

function fullsib(data::PopData; n::Int = 100, ploidy::Int = 2)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_parent(alleles, loc, ploidy = ploidy) for j in 1:2]
        [insertcols!(out_df, Symbol(prefix * "_off$j") => cross(p1, p2)) for j in 1:2]
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "fullsib")
    out_df.name = PooledArray(out_df.name)
    out_df.population = PooledArray(out_df.population)
    out_df.locus = PooledArray(out_df.locus)
    
    meta_df = DataFrame(
        :name => unique(out_df.name), 
        :population => "fullsib", 
        :ploidy => ploidy, 
        :longitude => Vector{Union{Missing, Float32}}(undef, n*2),
        :latitude => Vector{Union{Missing, Float32}}(undef, n*2),
        )
    return PopData(meta_df, out_df)
end


function halfsib(data::PopData; n::Int = 100, ploidy::Int = 2)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2,p3 = [simulate_parent(alleles, loc, ploidy = ploidy) for j in 1:3]
        insertcols!(out_df, Symbol(prefix * "_off1") => cross(p1, p2))
        insertcols!(out_df, Symbol(prefix * "_off2") => cross(p1, p3))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "halfsib")
    out_df.name = PooledArray(out_df.name)
    out_df.population = PooledArray(out_df.population)
    out_df.locus = PooledArray(out_df.locus)
    
    meta_df = DataFrame(
        :name => unique(out_df.name), 
        :population => "halfsib", 
        :ploidy => ploidy, 
        :longitude => Vector{Union{Missing, Float32}}(undef, n*2),
        :latitude => Vector{Union{Missing, Float32}}(undef, n*2),
        )
    return PopData(meta_df, out_df)
end


function unrelated(data::PopData; n::Int = 100, ploidy::Int = 2)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_parent(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_off1") => Tuple.(p1))
        insertcols!(out_df, Symbol(prefix * "_off2") => Tuple.(p2))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "unrelated")
    out_df.name = PooledArray(out_df.name)
    out_df.population = PooledArray(out_df.population)
    out_df.locus = PooledArray(out_df.locus)
    
    meta_df = DataFrame(
        :name => unique(out_df.name), 
        :population => "unrelated", 
        :ploidy => ploidy, 
        :longitude => Vector{Union{Missing, Float32}}(undef, n*2),
        :latitude => Vector{Union{Missing, Float32}}(undef, n*2),
        )
    return PopData(meta_df, out_df)
end