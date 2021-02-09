module  TestUtils

using PopGen
using PopGenSims
using Test

cats = @nancycats ;
cats2 = deepcopy(cats)

@testset "append PopData" begin
    @test length(append(cats, cats2).meta.name) == 2 * length(cats.meta.name)
    @test length(append(cats, cats2).loci.name) == 2 * length(cats.loci.name)
    @test typeof(append(cats, cats2)) == PopData
    append!(cats2, cats)
    @test length(cats2.meta.name) == 2 * length(cats.meta.name)
    @test length(cats2.loci.name) == 2 * length(cats.loci.name)
    @test typeof(cats2) == PopData 
end

@testset "allele pool" begin
    @test PopGenSims.allele_pool(cats.loci.genotype[1:30]) |> length == 56
    @test PopGenSims.allele_pool(cats.loci.genotype[1:30]) |> typeof <: NTuple
    a,b = PopGenSims.allele_pool(cats)
    @test eltype(a) == String
    @test length(a) == 9
    @test typeof(b) <: Dict
    @test length(b) == length(a)
    c = PopGenSims.simulate_sample(b, a, ploidy = 2)
    @test length(c) == length(a) == length(b)
    @test all(length.(c) .== 2)
    c3 = PopGenSims.simulate_sample(b, a, ploidy = 3)
    @test all(length.(c3) .== 3)
end

end # module