#author: Daniel Drapala 244939
#finding smallest x that x*(1/x)!=1
#x- number to start from (and next numbers are nextfloat())
function dividing(x)
    
    while Float64(x*(Float64(1.0)/x)) == Float64(1.0) && x < Float64(2.0)
        x = nextfloat(x)
    end
    return x
end
println("The smallest x which fulfill the equation is: $(dividing(Float64(nextfloat(1.0))))")
println("The other x which fulfill the equation is: $(dividing(Float64(nextfloat(1.5))))")


