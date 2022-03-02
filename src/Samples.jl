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
    Tuple.(sort.([sample(k, Weights(v), ploidy) for i in 1:n]))
end

"""
    simulate(data::PopData; n::Int)
Simulate `n` number of individuals per population using per-population
allele frequencies derived from a `PopData` object. Returns a new `PopData` object with `n` * `n_populations` samples.

    simulate(data::PopData; scale::Int)
Simulate individuals per population in the same proportions they appear in the PopData
using per-population allele frequencies. Simulation volume can be multiplied using `scale`,
i.e. if you want to keep the same proportions but generate twice the number of samples, `scale`
would be `2`. Returns a new `PopData` object with `n_samples` * `scale` samples.    

**Example**
```julia
julia> cats = @nanycats;

julia> sims = simulate(cats, n = 100)
PopData{Diploid, 9 Microsatellite Loci}
  Samples: 1700
  Populations: 17
  
julia> sims_prop = simulate(cats, scale = 3)
  PopData{Diploid, 9 Microsatellite Loci}
    Samples: 711
    Populations: 17
```
"""
function simulate(data::PopData; n::Int=0, scale::Int = 0)
    n == scale == 0 && throw(ArgumentError("Please use one of n (flat) or scale (proportional) keywords for simulations. See ?simulate for more info."))
    !iszero(n) & !iszero(scale) && throw(ArgumentError("Must use only one of n (flat) or scale (proportional) keywords for simulations. See ?simulate for more info."))
    !iszero(n) & iszero(scale) && return _simulateflat(data, n)
    iszero(n) & !iszero(scale) && return _simulatescale(data, scale)
end

function _simulateflat(data::PopData, n::Int)
    ploidy = data.metadata.ploidy
    length(ploidy) > 1 && error("Simulations do not work on mixed-ploidy data (yet)")
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

function _simulatescale(data::PopData, n::Int)
    ploidy = data.metadata.ploidy
    length(ploidy) > 1 && error("Simulations do not work on mixed-ploidy data (yet)")
    pops = unique(data.sampleinfo.population)
    popcounts =  [count(i -> i == j, data.sampleinfo.population) for j in pops]
    nloci = data.metadata.loci

    # instantiate output df
    simnames = repeat(["sim_" * "$i" for i in 1:(data.metadata.samples * n)], inner = nloci)
    popnames = reduce(vcat, [fill(i, nloci * j * n) for (i,j) in zip(pops,popcounts)])
    locinames = repeat(unique(data.genodata.locus), outer = (n * data.metadata.samples))
    geno_out = DataFrame(:name => simnames, :population => popnames, :locus => locinames, :genotype => similar(data.genodata.genotype, length(locinames)))
    # generate allele freqs per population
    gdf = groupby(data.genodata, [:population, :locus])
    freqs = DataFrames.combine(
        gdf,
        :genotype => allelefreq => :frq
    )
    freqs[:, :n] .= repeat(popcounts, inner = nloci)    
    transform!(freqs, [:frq, :n] => ((i,j) -> sample_locus.(i, j*n, ploidy)) => :frq)
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