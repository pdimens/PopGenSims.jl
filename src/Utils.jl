export append, append!

"""
    append!(data::PopData, data2::PopData)
Add the rows of `data2` to the end of `data`. This will add the samples present
in the second `PopData` object to the first `PopData` object (mutating it). 
**Note** that this is a simple appending, and you risk corrupting your `PopData` if
the two `PopData` objects do not have identical loci.

**Example**
```
julia> cats = @nancycats
PopData Object
  9 Microsatellite Loci
  Ploidy: 2
  Samples: 237
  Populations: 17
  Coordinates: absent

julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData Object
  9 Microsatellite Loci
  Ploidy: 2
  Samples: 100
  Populations: 1
  Coordinates: absent

julia> append!(cats, purrfect_pairs);

julia> cats
PopData Object
  9 Microsatellite
  Ploidy: 2
  Samples: 337
  Populations: 18
  Coordinates: absent
```
"""
function Base.append!(data::PopData, data2::PopData)
    if "parents" ∉ names(data.meta) && "parents" ∈ names(data2.meta)
        len = length(data.meta.name)
        insertcols!(
            data.meta, 
            :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
            )
    elseif "parents" ∉ names(data2.meta) && "parents" ∈ names(data.meta)
        len = length(data2.meta.name)
        insertcols!(
            data2.meta, 
            :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
            )
    end
    
    append!(data.meta, data2.meta)

    append!(data.loci, data2.loci)
    return data
end


"""
    append(data::PopData, data2::PopData)
Add the rows of `data2` to the end of `data`. This will combine the samples present
in both `PopData` objects and return a new `PopData` object. **Note** that this is 
a simple appending, and you risk corrupting your `PopData` if the two `PopData` 
objects do not have identical loci.

**Example**
```
julia> cats = @nanycats
PopData Object
  9 Microsatellite Loci
  Ploidy: 2
  Samples: 237
  Populations: 17
  Coordinates: absent

julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData Object
  9 Microsatellite Loci
  Ploidy: 2
  Samples: 100
  Populations: 1
  Coordinates: absent

julia> merged_cats = append(cats, purrfect_pairs)
PopData Object
  9 Microsatellite Loci
  Ploidy: 2
  Samples: 337
  Populations: 18
  Coordinates: absent
```
"""
function append(data::PopData, data2::PopData)
    tmp = PopData(copy(data.meta), copy(data.loci))
    tmp2 = data2
    if "parents" ∉ names(data.meta) && "parents" ∈ names(data2.meta)
        len = length(tmp.meta.name)
        insertcols!(
            tmp.meta, 
            :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
            )
    elseif "parents" ∉ names(data2.meta) && "parents" ∈ names(data.meta)
        tmp2 = PopData(copy(data2.meta), copy(data2.loci))
        len = length(tmp2.meta.name)
        insertcols!(
            tmp2.meta, 
            :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
            )
    end
    
    append!(tmp.meta, tmp2.meta)

    append!(tmp.loci, tmp2.loci)
    return tmp
end


function allele_pool(locus::T) where T <: GenoArray
  Tuple(Base.Iterators.flatten(skipmissing(locus)))
end

function allele_pool(data::PopData)
  # index dataframe by locus
  idx_df = groupby(data.loci, [:locus])
  # instantiate dict to store alleles
  #allele_dict = Dict{String,Tuple}()
  # pull out loci names
  loc = getindex.(keys(idx_df), :locus)
  allele_dict = Dict(i => allele_pool(idx_df[(;locus = i)].genotype) for i in loc)
  #[allele_dict[i] = allele_pool(idx_df[(;locus = i)].genotype) for i in loc]
  return string.(loc), allele_dict
end

"""
```
simulate_sample(alleles::Dict{String,NTuple}, loc::Vector{String}; ploidy::Int)
```
Using a global allele pool given by a Dict{loci,alleles} and a list of loci (`loc`), simulate
an individual with a given `ploidy`. Returns a Vector of genotypes.

**Example**
```
julia> cats = @nanycats ;
julia> loc, alleles = allele_pool(cats) ;
julia> simulate_sample(alleles, loc, ploidy = 2)
9-element Array{Array{Int16,1},1}:
 [139, 129]
 [146, 146]
 [145, 141]
 [126, 126]
 [150, 148]
 [148, 140]
 [185, 199]
 [91, 113]
 [208, 208]
```
"""
function simulate_sample(alleles::Dict{String,<:Tuple}, loc::Vector{String}; ploidy::Int)
  map(i -> rand(alleles[i], ploidy) ,loc)
end
