module  TestSibship

using PopGenCore
using PopGenSims
using Test

cats = @nancycats ;

@testset "Sibship Simulations" begin
    @testset "all sibship" begin
        sims_all = simulate_sibship(cats, fullsib = 50, halfsib = 50, unrelated = 50, parentoffspring = 50)
        @test sims_all isa PopData
        @test sims_all.metadata.samples == 400
        @test length(sims_all.genodata.name) == 3600
    end

    @testset "unrelated" begin
        sims_un = simulate_sibship(cats, unrelated = 50)
        @test sims_un isa PopData
        @test sims_un.metadata.samples == 100
        @test length(sims_un.genodata.name) == 900
    end

    @testset "halfsib" begin
        sims_hs = simulate_sibship(cats, halfsib = 50)
        @test sims_hs isa PopData
        @test sims_hs.metadata.samples == 100
        @test length(sims_hs.genodata.name) == 900
    end

    @testset "fullsib" begin
        sims_fs = simulate_sibship(cats, fullsib = 50)
        @test sims_fs isa PopData
        @test sims_fs.metadata.samples == 100
        @test length(sims_fs.genodata.name) == 900
    end

    @testset "parent offspring" begin
        sims_po = simulate_sibship(cats, parentoffspring = 50)
        @test sims_po isa PopData
        @test sims_po.metadata.samples == 100
        @test length(sims_po.genodata.name) == 900
    end
end

end # module