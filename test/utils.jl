module  TestUtils

using PopGenCore
using PopGenSims
using Test

cats = @nancycats ;
cats2 = copy(cats)

@testset "append PopData" begin
    tmp = append(cats, cats2)
    @test tmp.metadata.samples == 2 * cats.metadata.samples
    @test length(tmp.genodata.name) == 2 * length(cats.genodata.name)
    @test tmp isa PopData
    append!(cats2, cats)
    @test cats2.metadata.samples == 2 * cats.metadata.samples
    @test length(cats2.genodata.name) == 2 * length(cats.genodata.name)
    @test cats2 isa PopData 
end

@testset "allele pool" begin
    @test PopGenSims.allele_pool(cats.genodata.genotype[1:30]) |> length == 56
    @test PopGenSims.allele_pool(cats.genodata.genotype[1:30]) |> typeof <: NTuple
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