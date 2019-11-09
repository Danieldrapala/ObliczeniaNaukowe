#Daniel Drapla 244939

# T - type, p - starting point, r - given constant, n - number of iterations
function logical_iterate(T,p,r,n)
    P = zeros(T, n)
    np = T(p)
    for i = 1 : n
        np = p + r * p * (T(1.0) - p)
        p = np
        P[i] = p
    end
    return P
end

T = Float32
println("without mod")
A = logical_iterate(T,T(0.01),3,40)
for i = 1:length(A)
        println("$(A[i])")
    
end
println("with mod")
B = logical_iterate(T,T(0.01),3,10)
B[10] = round(B[10],digits=3,RoundDown)        #cuts to third place after comma
B1 = logical_iterate(T,B[10],3,40-10)
append!(B,B1)

for i = 1:length(B)
     
        println(" $(B[i])")
    
end

T = Float64
println("Float64")
C = logical_iterate(T,T(0.01),3,40)
for i = 1:length(C)
        println("$(C[i])")
    
end

