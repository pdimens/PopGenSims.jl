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
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 237
  Populations: 17

julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 100
  Populations: 1

julia> append!(cats, purrfect_pairs);

julia> cats
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 337
  Populations: 18
```
"""
function Base.append!(data::PopData, data2::PopData)
  n1 = data.metadata.samples
  pl1 = data.metadata.ploidy
  if "parents" ∉ names(data.sampleinfo) && "parents" ∈ names(data2.sampleinfo)
      len = data.metadata.samples
      insertcols!(
          data.sampleinfo, 
          :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
          )
  elseif "parents" ∉ names(data2.sampleinfo) && "parents" ∈ names(data.sampleinfo)
      len = length(data2.sampleinfo.name)
      insertcols!(
          data2.sampleinfo, 
          :parents => Vector{Union{Missing, Tuple{String,String}}}(undef, len)
          )
  end
  
  append!(data.sampleinfo, data2.sampleinfo)

  append!(data.genodata, data2.genodata)
  # update metadata
  data.metadata.samples = n1 + data2.metadata.samples
  if pl1 != data2.metadata.ploidy
    data.metadata.ploidy = Int8[pl1, data2.metadata.ploidy]
  end
  data.metadata.populations = length(unique(data.sampleinfo.population))
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
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 237
  Populations: 17


julia> purrfect_pairs = cross(cats, "N200", "N7", generation = "F1")
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 100
  Populations: 1

julia> merged_cats = append(cats, purrfect_pairs)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 337
  Populations: 18
```
"""
function append(data::PopData, data2::PopData)
  tmp  = copy(data)
  append!(tmp, data2)
  return tmp
end


function allele_pool(locus::T) where T <: GenoArray
  Tuple(Base.Iterators.flatten(skipmissing(locus)))
end

function allele_pool(data::PopData)
  # index dataframe by locus
  idx_df = groupby(data.genodata, [:locus])
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
function simulate_sample(alleles::Dict{String,<:Tuple}, loc::Vector{String}; ploidy::Signed)
  map(i -> rand(alleles[i], ploidy) ,loc)
end
