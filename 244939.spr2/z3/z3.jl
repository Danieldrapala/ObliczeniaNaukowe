# Daniel Drapała 244939 
# 4.11.2019


include("hilb.jl")
include("matcond.jl")

# A - matrix, b - vector machx - calculated vector, x - real vector

gaussian_elimination(A,b) = A\b
inversion(A,b) = inv(A)*b
approx_error(machx, x) = norm(machx - x) / norm(x)

function print_result(T,A,n)
    x = ones(T, n)
    b = A*x
    gauss_x = gaussian_elimination(A,b)
    inv_x = inversion(A,b)
    gauss_err = approximation_error(gauss_x,x)
    inv_err = approximation_error(inv_x,x)

    # println("$(n)x$(n) & $(rank(A)) & $(cond(A)) & $(gauss_err) & $(inv_err) \\\\")
    # @printf("%dx%d & %d & %.15e & %.15e & %.15e \\\\\n",n,n,rank(A),cond(A),gauss_err,inv_err)
    println("Size: $(n)x$(n)\tRank: $(rank(A))")
    println("Cond: $(cond(A))")
    println("Gauss error: $(gauss_err)\nInversion error: $(inv_err)\n")
end

T = Float64
# maximum size for Hilbert's matrix
matrix_max = 20
# conds for random matrix
c = [T(1), T(10), T(10^3), T(10^7), T(10^12), T(10^16)]
# sizes for random matrix
n = [5, 10, 20]

# Running program for Hilbert's matrices of different sizes
println("Hilbert Matrix\n")
for i = 1: matrix_max
    print_result(T,hilb(i),i)
end

# Running program for random matrices of different size and cond
println("\n\nRandom Matrix\n")
for i = 1: length(n)
    for j = 1: length(c)
        print_result(T, matcond(n[i],c[j]), n[i])
    end
end