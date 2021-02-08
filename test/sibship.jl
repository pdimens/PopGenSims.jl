module  TestSibship

using PopGen
using PopGenSims
using Test

cats = @nancycats ;

@testset "unrelated" begin
    sims = simulate_sibship(cats, n = 50, relationship= "unrelated")
    @test typeof(sims) == PopData
    @test length(sims.meta.name) == 100
    @test length(sims.loci.name) == 900
end

@testset "halfsib" begin
    sims = simulate_sibship(cats, n = 50, relationship= "halfsib")
    @test typeof(sims) == PopData
    @test length(sims.meta.name) == 100
    @test length(sims.loci.name) == 900
end

@testset "fullsib" begin
    sims = simulate_sibship(cats, n = 50, relationship= "fullsib")
    @test typeof(sims) == PopData
    @test length(sims.meta.name) == 100
    @test length(sims.loci.name) == 900
end

@testset "parent-offspring" begin
    sims = simulate_sibship(cats, n = 50, relationship= "parent-offspring")
    @test typeof(sims) == PopData
    @test length(sims.meta.name) == 100
    @test length(sims.loci.name) == 900
end


end # module