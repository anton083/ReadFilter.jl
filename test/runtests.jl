using ReadFilter
using Test

@testset "ReadFilter.jl" begin
    @testset "io.jl" begin
        @test sequence_subranges(10, 5, 2) == [1:5, 4:8, 6:10]
    end
end