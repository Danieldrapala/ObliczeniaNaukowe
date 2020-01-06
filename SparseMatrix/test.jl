#Daniel Drapala 244939 - author
  

include("blocksys.jl")
include("fileinout.jl")

using Test
using Main.blocksys
using Main.fileinout

# Default test sizes
sizes = [16, 10000, 50000]

@testset "Features test for matrix $size" for size in sizes
    A, n, l = load_matrix_ffile("./Dane$(size)_1_1/A.txt")
    b = load_vector_ffile("./Dane$(size)_1_1/b.txt")
    
    @testset "Computing b from A" begin
        @test isapprox(right_side_vector(A, n, l), b)
    end
    A, n, l = load_matrix_ffile("./Dane$(size)_1_1/A.txt")
    b = load_vector_ffile("./Dane$(size)_1_1/b.txt")

    @testset "Gaussian elimination" begin
        @test isapprox(\(A, b),solving_eq_after_GE(A, n, l, b))
    end
    A, n, l = load_matrix_ffile("./Dane$(size)_1_1/A.txt")
    b = load_vector_ffile("./Dane$(size)_1_1/b.txt")

    @testset "Pivoted gaussian elimination" begin
        @test isapprox( \(A, b),solving_eq_after_GEWP(A, n, l, b))
    end
    A, n, l = load_matrix_ffile("./Dane$(size)_1_1/A.txt")
    b = load_vector_ffile("./Dane$(size)_1_1/b.txt")

    @testset "LU" begin
        @test isapprox(\(A, b),solving_eq_after_LU(A, n, l, b))
    end
    A, n, l = load_matrix_ffile("./Dane$(size)_1_1/A.txt")
    b = load_vector_ffile("./Dane$(size)_1_1/b.txt")
    @testset "Pivoted LU" begin
        @test isapprox(\(A, b),solving_eq_after_LUWP(A, n, l, b) )
    end
end