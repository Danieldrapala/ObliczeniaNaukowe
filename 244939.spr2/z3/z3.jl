# Daniel Drapała 244939 
# 4.11.2019


include("hilb.jl")
include("matcond.jl")

# A - matrix, b - vector, x - real vector

function print_values(T,A,n)
    #create vector x with ones
    x = ones(T, n)
    b = A*x
    gauss_x = A\b
    inv_x = inv(A)*b
    gauss_err = norm(gauss_x - x) / norm(x)
    inv_err = norm(inv_x - x) / norm(x)
    print(" Size:= $(n)x$(n) = Rank = $(rank(A)) =Cond: = $(cond(A)) =Gauss error: = $(gauss_err) =Inversion error:=  $(inv_err) =\n")
end

T = Float64

println("Hilbert Matrix")
for i = 3: 20
    print_values(T,hilb(i),i)
end

c = [T(1), T(10), T(10^3), T(10^7), T(10^12), T(10^16)]
n = [5, 10, 20]

println("Random Matrix")
for i = 1: length(n)
    for j = 1: length(c)
        print_values(T, matcond(n[i],c[j]), n[i])
    end
end

