
#finding smallest x that x*(1/x)!=1
x = nextfloat(Float64(1.0))

while Float64(x*(Float64(1.0)/x))==Float64(1.0) && x<Float64(2.0)
    global x
    x = nextfloat(x)
end
println( "The smallest x is $(x) ")
