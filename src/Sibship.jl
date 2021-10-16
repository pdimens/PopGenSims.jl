export simulate_sibship

"""
    _cross(parent1::Vector{Vector{T}}, parent2::Vector{Vector{T}}) where T <: Signed
Simulate a mating cross between two parents, generating one offspring with the same
ploidy as `parent1`. This variant of `cross` is used internally for `simulate_sibship`.
"""
function _cross(parent1::Vector{Vector{T}}, parent2::Vector{Vector{T}}) where T <: Signed
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


function _parentoffspring(alleles::Dict, loc::Vector{String}, n::Int, ploidy::Signed, padding::Int)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim" * lpad(i, padding, '0')
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_parent") =>  Tuple.(p1))
        insertcols!(out_df, Symbol(prefix * "_offspring") => _cross(p1, p2))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "parent_offspring")
    #out_df.name = PooledArray(out_df.name, compress = true)
    #out_df.population = PooledArray(out_df.population, compress = true)
    #out_df.locus = PooledArray(out_df.locus, compress = true) 
    return out_df
end


function _fullsib(alleles::Dict, loc::Vector{String}, n::Int, ploidy::Signed, padding::Int)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim" * lpad(i, padding, '0')
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        [insertcols!(out_df, Symbol(prefix * "_fullsib_$j") => _cross(p1, p2)) for j in 1:2]
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "fullsib")
    #out_df.name = PooledArray(out_df.name, compress = true)
    #out_df.population = PooledArray(out_df.population, compress = true)
    #out_df.locus = PooledArray(out_df.locus, compress = true)
    return out_df
end



function _halfsib(alleles::Dict, loc::Vector{String}, n::Int, ploidy::Signed, padding::Int)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim" * lpad(i, padding, '0')
        p1,p2,p3 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:3]
        insertcols!(out_df, Symbol(prefix * "_halfsib_1") => _cross(p1, p2))
        insertcols!(out_df, Symbol(prefix * "_halfsib_2") => _cross(p1, p3))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "halfsib")
    #out_df.name = PooledArray(out_df.name, compress = true)
    #out_df.population = PooledArray(out_df.population, compress = true)
    #out_df.locus = PooledArray(out_df.locus, compress = true)
    return out_df
end


function _unrelated(alleles::Dict, loc::Vector{String}, n::Int, ploidy::Signed, padding::Int)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim" * lpad(i, padding, '0')
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_unrelated_1") => Tuple.(p1))
        insertcols!(out_df, Symbol(prefix * "_unrelated_2") => Tuple.(p2))
    end
    out_df = rename!(select!(DataFrames.stack(out_df, Not(:locus)), 2, 1, 3), [:name, :locus, :genotype])
    insertcols!(out_df, 2, :population => "unrelated")
    #out_df.name = PooledArrayout_df.name, compress = true)
    #out_df.population = PooledArray(out_df.population, compress = true)
    #out_df.locus = PooledArray(out_df.locus, compress = true)
    return out_df
end

"""
    simulate_sibship(data::PopData; fullsib::Int, halfsib::Int, unrelated::Int, parentoffspring::Int, ploidy::Signed)
Simulate mating crosses to generate sample pairs with any combination of the specified relationships, 
returning a `PopData` object. The simulations will first generate parents of a given
`ploidy` (inferred or specified) by drawing alleles from a global allele pool derived
from the given `data` (i.e. weighted by their frequencies).

#### Relationship
Simulated parents will be crossed to generate offspring depending on the relationship:
- `fullsib` : 2 parents generate 2 full-sibling offspring, return 2 offspring
- `halfsib` : 3 parents generate 2 half-sibling offspring, returns 2 offspring
- `unrelated` : returns 2 randomly generated individuals from the global allele pools
- `parentoffspring` : 2 parents generate 1 offspring, returns 1 offspring and 1 parent

#### Identifying pairs
The relationship between the newly generated samples can be identified by:
- Sample `name`s will specify their simulation number, relationship, and whether parent or offspring
    - Naming convention: [simulation #]_[relationship]_[offspring #]
    - example: sim005_fullsib_1 = [simulation 005]_[full sibling]_[offspring 1]
- Their `population` name will be that of their relationship (e.g. "fullsib")

#### Ploidy
If the samples in your `PopData` are of a single ploidy, then `ploidy = 0` (the default) will infer the ploidy
and generate parents and offspring according to the ploidy of your data. If you have mixed-ploidy data or wish 
to generate parents and offspring of a ploidy different than the source `PopData` you can specify the ploidy
with which to simulate parents and offspring. For example, if your `PopData` is diploid, but you wish to generate
triploid or octoploid parents and offspring, you would specify `ploidy = 3` or `ploidy = 8` repectively. 
#### Odd ploidy
If trying to create offspring with an odd ploidy (3,5, etc.), each parent has a 50% chance of 
contributing (½ × ploidy) + 1 alleles for all loci to the offspring. In other words, if ploidy = 3,
there's a 50% chance parent_1 will give 2 alleles for every locus for that simulated cross.

**Example**
```
julia> cats = @nanycats ;

julia> cat_sims = simulate_sibship(cats, fullsib = 10, halfsib = 50)
PopData{Diploid, 9 Microsatellite loci}
  Samples: 120
  Populations: 2

julia> cat_sims.sampleinfo
120×3 DataFrame
 Row │ name             population  ploidy 
     │ String           String      Int64  
─────┼─────────────────────────────────────
   1 │ sim01_fullsib_1  fullsib          2
   2 │ sim01_fullsib_2  fullsib          2
   3 │ sim02_fullsib_1  fullsib          2
   4 │ sim02_fullsib_2  fullsib          2
   5 │ sim03_fullsib_1  fullsib          2
   6 │ sim03_fullsib_2  fullsib          2
  ⋮  │        ⋮             ⋮         ⋮
 115 │ sim48_halfsib_1  halfsib          2
 116 │ sim48_halfsib_2  halfsib          2
 117 │ sim49_halfsib_1  halfsib          2
 118 │ sim49_halfsib_2  halfsib          2
 119 │ sim50_halfsib_1  halfsib          2
 120 │ sim50_halfsib_2  halfsib          2
                           108 rows omitted
```
"""
function simulate_sibship(data::PopData; fullsib::Int = 0, halfsib::Int = 0, unrelated::Int = 0, parentoffspring::Int = 0, ploidy::Signed = 0)
    if iszero(sum([fullsib, halfsib, unrelated, parentoffspring]))
        throw(ArgumentError("Please specify at least one of: \n- \"fullsib\" \n- \"halfsib\" \n- \"unrelated\"\n- \"parentoffspring\""))
    end
    # automatic ploidy finding
    if ploidy == 0
        ploids = data.metadata.ploidy
        if ploids isa AbstractVector
            error("For mixed ploidy data, please specify a single ploidy with which to generate parents and offspring")
        else
            ploidy += ploids    
        end
    end
    loc, alleles = allele_pool(data)
    # how many digits to pad the offspring names
    padding = length(string(maximum([fullsib, halfsib, unrelated, parentoffspring])))
    # perform the simulation if the integer > 0, otherwise return an empty boolean vector
    # the empty vector is just to skip over with Base.Iterators.filter
    fs = fullsib > 0 ? _fullsib(alleles, loc, fullsib, ploidy, padding) : Bool[]
    hs = halfsib > 0 ? _halfsib(alleles, loc, halfsib, ploidy, padding) : Bool[]
    unrl = unrelated > 0 ? _unrelated(alleles, loc, unrelated, ploidy, padding) : Bool[]
    poff = parentoffspring > 0 ? _parentoffspring(alleles, loc, parentoffspring, ploidy, padding) : Bool[]
    # combine the results together into a single df
    geno_df = reduce(vcat, Base.Iterators.filter(!isempty, (fs, hs, unrl, poff)))
    geno_df.name = PooledArray(geno_df.name, compress = true)
    geno_df.population = PooledArray(geno_df.population, compress = true)
    geno_df.locus = PooledArray(geno_df.locus, compress = true)
    #meta_df = select(unique(geno_df, :name), 1, 2)
    #insertcols!(meta_df, :ploidy => ploidy)
    return PopData(geno_df)
end