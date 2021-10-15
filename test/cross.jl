module  TestCross

using PopGenCore
using PopGenSims
using Test

cats = @nancycats ;

@testset "Crosses" begin
    @testset "sample_genotypes" begin
        @test length(PopGenSims.sample_genotype(cats.loci.genotype[30], 2)) == 2
        @test eltype(PopGenSims.sample_genotype(cats.loci.genotype[30],1)) == eltype(cats.loci.genotype[30])
        @test length(PopGenSims.sample_genotype(missing, 2)) == 1
        @test first(PopGenSims.sample_genotype(missing, 2)) === missing
    end

    @testset "crosses" begin
        f1 = cross(cats, "N111", "N107", n = 10)
        @test f1 isa PopData
        @test length(f1.meta.name) == 10
        @test length(f1.loci.name) == 90
        f1 = cross(cats, "N111", "N107", n = 10, generation = "firstgen")
        @test f1 isa PopData
        @test length(f1.meta.name) == 10
        @test length(f1.loci.name) == 90
        f2 = cross(cats => "N111", f1 => "firstgen_offspring_10", n = 10)
        @test f2 isa PopData
        @test length(f2.meta.name) == 10
        @test length(f2.loci.name) == 90
        f2 = cross(cats => "N111", f1 => "firstgen_offspring_10", n = 10, generation = "F2")
        @test f2 isa PopData
        @test length(f2.meta.name) == 10
        @test length(f2.loci.name) == 90
    end
end
end # module TestCross