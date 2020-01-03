
  
include("blocksys.jl")
include("fileinout.jl")

using Main.blocksys,Main.fileinout

# (A,n,l) = load_matrix_ffile("Dane10000_1_1/A.txt")
# b = load_vector_ffile("Dane10000_1_1/b.txt")
# # println(gaussian_elimination_with_pivots(A, n, l, b))
# p = matrix_to_LU_with_pivots(A, n, l)
# save_values("./wynikLUPtxt",solve_from_LU_with_pivots(A, n, l, b, p),n,true)

(A,n,l) = load_matrix_ffile("Dane10000_1_1/A.txt")
b = load_vector_ffile("Dane10000_1_1/b.txt")
matrix_to_LU(A,n,l)
save_values("./wynikLU.txt",solve_from_LU(A, n, l, b),n,true)

(A,n,l) = load_matrix_ffile("Dane10000_1_1/A.txt")
b = load_vector_ffile("Dane10000_1_1/b.txt")
save_values("./wynikGUP.txt",gaussian_elimination_with_pivots(A, n, l, b),n,true)

(A,n,l) = load_matrix_ffile("Dane10000_1_1/A.txt")
b = load_vector_ffile("Dane10000_1_1/b.txt")
save_values("./wynikGU.txt",gaussian_elimination(A, n, l, b),n,true)

(A,n,l) = load_matrix_ffile("Dane10000_1_1/A.txt")
b = right_side_vector(A,n,l)
save_values("./wynikGU1.txt",gaussian_elimination(A, n, l, b),n,true)