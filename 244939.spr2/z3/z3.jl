# Daniel Drapała 244939 
# 4.11.2019


include("hilb.jl")
include("matcond.jl")

# A - matrix, b - vector machx - calculated vector, x - real vector
gauss(A,b) = A\b
inver(A,b) = inv(A)*b
approx_err(machx, x) = norm(machx - x) / norm(x)

function print_result(T,A,n)
    #create vector x with ones
    x = ones(T, n)
    b = A*x
    gauss_x = gauss(A,b)
    inv_x = inver(A,b)
    gauss_err = approx_err(gauss_x,x)
    inv_err = approx_err(inv_x,x)
    print(" Size:= $(n)x$(n) = Rank = $(rank(A))")
    print(" =Cond: = $(cond(A)) =")
    print(" Gauss error: = $(gauss_err) =Inversion error:=  $(inv_err) =\n")
end

T = Float64

println("Hilbert Matrix")
for i = 3: 20
    print_result(T,hilb(i),i)
end

c = [T(1), T(10), T(10^3), T(10^7), T(10^12), T(10^16)]
n = [5, 10, 20]

println("Random Matrix")
for i = 1: length(n)
    for j = 1: length(c)
        print_result(T, matcond(n[i],c[j]), n[i])
    end
end