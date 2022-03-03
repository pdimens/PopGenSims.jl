module SamplesTest

using PopGenCore
using PopGenSims
using Test

cats = @nancycats ;
d = Dict("1" => 5, "4" => 7, "11" => 3)

@testset "Simulate Samples" begin
    @testset "sample_locus" begin
        d = Dict(133 => 0.125, 135 => 0.5625, 143 => 0.25, 137 => 0.0625)
        @test PopGenSims.sample_locus(d, 3, 2) |> length == 3
        @test PopGenSims.sample_locus(d, 3, 2) |> first |> length == 2
        @test PopGenSims.sample_locus(d, 3, 3) |> length == 3
        @test PopGenSims.sample_locus(d, 3, 3) |> first |> length == 3
    end

    @testset "simulate errors" begin
        @test_throws ArgumentError simulate(cats)
        @test_throws ArgumentError simulate(cats, n = 1, scale = 1)
        @test_throws ArgumentError simulate(cats, n = d, scale = 2)
    end

    @testset "simulate flat" begin
        sims = simulate(cats , n = 100)
        @test sims isa PopData
        @test sims.metadata.samples == 1700
        @test length(sims.genodata.name) == 15300
    end

    @testset "simulate proportional" begin
        sims = simulate(cats , scale = 1)
        @test sims isa PopData
        @test sims.metadata.populations == cats.metadata.populations
        @test sims.metadata.samples == cats.metadata.samples
        @test length(sims.genodata.name) == length(cats.genodata.name)
        sims = simulate(cats , scale = 4)
        @test sims isa PopData
        @test sims.metadata.populations == cats.metadata.populations
        @test sims.metadata.samples == cats.metadata.samples * 4
        @test length(sims.genodata.name) == length(cats.genodata.name) * 4
    end

    @testset "simulate arbitrary" begin
        sims = simulate(cats , n = d)
        @test sims isa PopData
        @test sims.metadata.samples == sum(values(d))
        @test length(sims.genodata.name) == sum(values(d)) * cats.metadata.loci
        @test length(d) == sims.metadata.populations
    end
end

end # module