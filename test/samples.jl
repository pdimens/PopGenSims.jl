module SamplesTest

using PopGenCore
using PopGenSims
using Test

cats = @nancycats ;

@testset "Simulate Samples" begin
    @testset "sample_locus" begin
        d = Dict(133 => 0.125, 135 => 0.5625, 143 => 0.25, 137 => 0.0625)
        @test PopGenSims.sample_locus(d, 3, 2) |> length == 3
        @test PopGenSims.sample_locus(d, 3, 2) |> first |> length == 2
        @test PopGenSims.sample_locus(d, 3, 3) |> length == 3
        @test PopGenSims.sample_locus(d, 3, 3) |> first |> length == 3
    end



    @testset "simulate" begin
        sims = simulate(cats , n = 100)
        @test sims isa PopData
        @test length(sims.meta.name) == 1700
        @test length(sims.loci.name) == 15300
    end
end

end # module