#author Daniel Drapala 244939
#delta- distance between two nextfloat numbers
#from - function starts to get nextfloat from this number
#steps - it describes how many next float numbers function has to find 
function stepforward(delta, from, steps)
    fl = Float64(from)
    for i = 1 : steps
        fl += delta
        println("$(bitstring(fl))")
        
    end
end
function stepbackward(delta, from, steps)
    fl = Float64(from)
    for i = 1 : steps
        fl -= delta
        
        println("$(bitstring(fl))")
        
    end
end


stepforward(Float64(2^-53),0.5,4)
stepbackward(Float64(2^-53),1,4)

stepforward(Float64(2^-51),2.0,4)
stepbackward(Float64(2^-51),4.0,4)

