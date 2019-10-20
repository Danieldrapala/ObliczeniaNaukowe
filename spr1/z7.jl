#author: Daniel Drapala 244939

#f - function
#x - argument of function
#h - very small factor

f(x) = sin(x) + cos(3x)
der(x) = cos(x) - 3sin(3x)

approx_der(f, x, h) = (f(x+h)-f(x))/h
approx_err(val, approx) = abs(val - approx)

for i = 0 : 54
    h = Float64(2.0)^(-i)
    rv = der(Float64(1.0))
    av = approx_der(f, Float64(1.0), h)
    err = approx_err(rv, av)
    println("2-$(i)      $(rv)  $(av)     $(err)      $(1+h)    ")
end
