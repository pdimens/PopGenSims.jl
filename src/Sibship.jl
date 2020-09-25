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


function parentoffspring(data::PopData; n::Int, ploidy::Int)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_parent") =>  Tuple.(p1))
        insertcols!(out_df, Symbol(prefix * "_offspring") => _cross(p1, p2))
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


function fullsib(data::PopData; n::Int, ploidy::Int)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        [insertcols!(out_df, Symbol(prefix * "_fs_off$j") => _cross(p1, p2)) for j in 1:2]
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


function halfsib(data::PopData; n::Int, ploidy::Int)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2,p3 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:3]
        insertcols!(out_df, Symbol(prefix * "_hs_off1") => _cross(p1, p2))
        insertcols!(out_df, Symbol(prefix * "_hs_off2") => _cross(p1, p3))
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


function unrelated(data::PopData; n::Int, ploidy::Int)
    loc, alleles = allele_pool(data)
    out_df = DataFrame(:locus => loc)
    for i in 1:n
        prefix = "sim$i"
        p1,p2 = [simulate_sample(alleles, loc, ploidy = ploidy) for j in 1:2]
        insertcols!(out_df, Symbol(prefix * "_un_off1") => Tuple.(p1))
        insertcols!(out_df, Symbol(prefix * "_un_off2") => Tuple.(p2))
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

"""
    simulate_sibship(data::PopData; n::Int, relationship::String, ploidy::Int)
Simulate mating crosses to generate `n` sample pairs (default: `500`) having the specified `relationship`, 
returning a `PopData` object. The simulations will first generate parents of a given `ploidy` (inferred or specified) 
by drawing alleles from a global allele pool derived from the given `data` (i.e. weighted by their frequencies).

#### Relationship
Simulated parents will be crossed to generate offspring depending on the `relationship`:
- `"fullsib"` : 2 parents generate 2 full-sibling offspring, return 2 offspring
- `"halfsib` : 3 parents generate 2 half-sibling offspring, returns 2 offspring
- `"unrelated"` : returns 2 randomly generated individuals from the global allele pools
- `"parent-offspring"` : 2 parents generate 1 offspring, returns 1 offspring and 1 parent

#### Identifying pairs
The relationship between the newly generated samples can be identified by:
- Sample `name`s will specify their simulation number, relationship, and whether parent or offspring
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
julia> cats = nancycats() ;

julia> fullsib_sims = simulate_sibship(cats, n = 50, relationship= "fullsib")
PopData Object
  Markers: Microsatellite
  Ploidy: 2
  Samples: 100
  Loci: 9
  Populations: 1
  Longitude: absent
  Latitude: absent

julia> fullsib_sims.meta_df100×5 DataFrame
│ Row │ name          │ population │ ploidy │ longitude │ latitude │
│     │ String        │ String     │ Int64  │ Float32?  │ Float32? │
├─────┼───────────────┼────────────┼────────┼───────────┼──────────┤
│ 1   │ sim1_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
│ 2   │ sim1_fs_off2  │ fullsib    │ 2      │ missing   │ missing  │
│ 3   │ sim2_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
│ 4   │ sim2_fs_off2  │ fullsib    │ 2      │ missing   │ missing  │
│ 5   │ sim3_fs_off1  │ fullsib    │ 2      │ missing   │ missing  │
⋮
│ 95  │ sim48_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 96  │ sim48_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
│ 97  │ sim49_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 98  │ sim49_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
│ 99  │ sim50_fs_off1 │ fullsib    │ 2      │ missing   │ missing  │
│ 100 │ sim50_fs_off2 │ fullsib    │ 2      │ missing   │ missing  │
```
"""
function simulate_sibship(data::PopData; n::Int = 500, relationship::String = "nothing", ploidy::Int = 0)
    if relationship == "nothing"
        error("Please use the keyword \'relationship\' and specify one of: \n- \"fullsib\" \n- \"halfsib\" \n- \"unrelated\"\n- \"parent-offspring\"") 
    elseif relationship ∉ ["fullsib", "halfsib", "unrelated", "parent-offspring"]
        error("relationship = \"$relationship\" is invalid, please specify one of: \n- \"fullsib\" \n- \"halfsib\" \n- \"unrelated\"\n- \"parent-offspring\"")
    end
    # automatic ploidy finding
    if ploidy == 0
        ploids = unique(data.meta.ploidy)
        if length(ploids) != 1
            error("For mixed ploidy data, please specify a single ploidy with which to generate parents and offspring")
        else
            ploidy += first(ploids)    
        end
    end
    if relationship == "fullsib"
        fullsib(data, n = n, ploidy = ploidy)
    elseif relationship == "halfsib"
        halfsib(data, n = n, ploidy = ploidy)
    elseif relationship == "unrelated"
        unrelated(data, n = n, ploidy = ploidy)
    else relationship == "halfsib"
        parentoffspring(data, n = n, ploidy = ploidy)
    end
end