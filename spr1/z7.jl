#author: Daniel Drapala 244939

#f - function
#x - argument of function
#h - very small factor

f(x) = sin(x) + cos(3x)
deltaf(x) = cos(x) - 3sin(3x)

approximated_derivative(f, x, h) = (f(x+h)-f(x))/h
approximation_error(real_val, approximated_val) = abs(real_val - approximated_val)

for i = 0 : 54
    h = Float64(2.0)^(-i)
    rv = deltaf(Float64(1.0))
    av = approximated_derivative(f, Float64(1.0), h)
    err = approximation_error(rv, av)
    println("h = 2^-$(i)\tAppr: $(av)\tReal: $(rv)\tErr: $(err)\t 1+h = $(1+h)\n")
end
