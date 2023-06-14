using ReadFilter
using Test

@testset "ReadFilter.jl" begin
    @testset "io.jl" begin
        @test get_subranges(10, 5, 2) == [1:5, 4:8, 6:10]
        @test get_subranges(10, 10, 2) == [1:10]
    end
end