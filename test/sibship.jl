module  TestSibship

using PopGen
using PopGenSims
using Test

cats = @nancycats ;
sims_un = simulate_sibship(cats, n = 50, relationship= "unrelated")
sims_hs = simulate_sibship(cats, n = 50, relationship= "halfsib")
sims_fs = simulate_sibship(cats, n = 50, relationship= "fullsib")
sims_po = simulate_sibship(cats, n = 50, relationship= "parent-offspring")


@testset "unrelated" begin
    @test typeof(sims_un) == PopData
    @test length(sims_un.meta.name) == 100
    @test length(sims_un.loci.name) == 900
end

@testset "halfsib" begin
    @test typeof(sims_hs) == PopData
    @test length(sims_hs.meta.name) == 100
    @test length(sims_hs.loci.name) == 900
end

@testset "fullsib" begin
    @test typeof(sims_fs) == PopData
    @test length(sims_fs.meta.name) == 100
    @test length(sims_fs.loci.name) == 900
end

@testset "parent offspring" begin
    @test typeof(sims_po) == PopData
    @test length(sims_po.meta.name) == 100
    @test length(sims_po.loci.name) == 900
end


end # module