#author: Daniel Drapala 244939

#x-argument
#f-function f
#g-function g
    for i = 1 : 9
        x = Float64(8)^(-i)
        f = sqrt(x^2+1)-1
        g = x^2/(sqrt(x^2+1)+1)
        println("$(g)")
    end

