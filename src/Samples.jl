export simulate

"""
    sample_locus(locus::Dict, n::Int, ploidy::Signed)
Internal function used by `simulate` to take a `Dict` of alleles => frequencies of a locus and return
`n` number of genotypes (n_alleles = `ploidy`) by using weighted sampling of the
allele-frequency pairs.

**Example**
```julia
d = Dict(133 => 0.125,135 => 0.5625,143 => 0.25,137 => 0.0625)

julia> sample_locus(d, 3, 2)
5-element Array{Tuple{Int16,Int16},1}:
 (133, 135)
 (135, 135)
 (143, 137)

julia> sample_locus(d, 3, 3)
5-element Array{Tuple{Int16,Int16,Int16},1}:
 (135, 135, 133)
 (143, 135, 133)
 (137, 135, 135)
```
"""
function sample_locus(locus::Dict, n::Int, ploidy::Signed)
    isempty(locus) && return fill(missing, n)
    k,v = collect(keys(locus)), collect(values(locus))
    alleles = [sample(k, Weights(v), n) for i in 1:ploidy]
    Tuple.(sort.(eachrow(hcat(alleles...))))
end


"""
    simulate(data::PopData; n::Int = 100)
Simulate `n` number of individuals per population (default: `100`) using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object.

**Example**
```julia
cats = @nanycats;

julia> sims = simulate(cats , n = 100)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 1700
  Populations: 17
  
julia> sims.sampleinfo

  1700×5 DataFrame
  Row │ name      population  ploidy   
      │ String    String      Int8      
──────┼───────────────────────────────
    1 │ sim_1     1                2    
    2 │ sim_2     1                2    
    3 │ sim_3     1                2    
    4 │ sim_4     1                2    
    5 │ sim_5     1                2    
  ⋮   │    ⋮          ⋮         ⋮ 
 1697 │ sim_1697  17               2  
 1698 │ sim_1698  17               2  
 1699 │ sim_1699  17               2  
 1700 │ sim_1700  17               2  
                                         1691 rows omitted

julia> sims.genodata
15300×4 DataFrame
   Row │ name      population  locus   genotype   
       │ String    String      String  Tuple…?    
───────┼──────────────────────────────────────────
     1 │ sim_1     1           fca8    (135, 143)
     2 │ sim_1     1           fca23   (136, 146)
     3 │ sim_1     1           fca43   (141, 145)
     4 │ sim_1     1           fca45   (120, 126)
     5 │ sim_1     1           fca77   (156, 156)
   ⋮   │    ⋮          ⋮         ⋮         ⋮
 15297 │ sim_1700  17          fca78   (150, 150)
 15298 │ sim_1700  17          fca90   (197, 197)
 15299 │ sim_1700  17          fca96   (113, 113)
 15300 │ sim_1700  17          fca37   (208, 208)
                                15291 rows omitted
```
"""
function simulate(data::PopData; n::Int = 100)
    data.metadata.ploidy != 2 && error("Simulations do not work on mixed-ploidy data (yet)")
    ploidy = data.metadata.ploidy
    pops = unique(data.sampleinfo.population)
    npops = data.metadata.populations
    nloci = data.metadata.loci

    # instantiate output df
    simnames = repeat(["sim_" * "$i" for i in 1:(n*npops)], inner = nloci)
    popnames = repeat(pops, inner = (nloci * n))
    locinames = repeat(unique(data.genodata.locus), outer = (n * npops))
    geno_out = DataFrame(:name => simnames, :population => popnames, :locus => locinames, :genotype => similar(data.genodata.genotype, length(locinames)))

    # generate allele freqs per population
    gdf = groupby(data.genodata, [:population, :locus])
    freqs = DataFrames.combine(
        gdf,
        :genotype => allelefreq => :frq
    )
    # create new genotypes
    transform!(freqs, :frq => (i -> sample_locus.(i,n,ploidy)) => :frq)
    
    # populate out df
    out_gdf = groupby(geno_out, :population)
    geno_gdf = groupby(freqs, :population)
    for pop in pops
        out_gdf[(population = pop,)][:,:genotype] = reduce(hcat, geno_gdf[(population = pop,)].frq) |> permutedims |> vec
    end
    transform!(
        geno_out,
        :name => (i -> PooledArray(i, compress = true)) => :name,
        :population => (i -> PooledArray(i, compress = true)) => :population,
        :locus => (i -> PooledArray(i, compress = true)) => :locus,
        :genotype
    )
    PopData(geno_out)
end