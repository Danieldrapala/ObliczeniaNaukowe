#author Daniel Drapala 244939

#t- type of numbers
#x- vector n0 one
#y- vector n0 two
#temp- temporary array used for sorting numbers in vector
#every sum- is a partial sum like sum1 and sum2 have to be added to make scalar product
function algoc(t,x,y)
    temp  =t[]
    for i=1 : length(x)
        push!(temp,x[i]*y[i])
    end

    sort!(temp, rev = true)   #sorting array in reverse order (biggest to smallest)
    sum1 ::t=t(0)
    sum2 ::t=t(0)
    j=length(temp)
    for i = 1:length(temp)
        if temp[i] > 0
            sum1 += temp[i]  
        end
        if temp[j]<0
            sum2+=temp[j]
        end
        j-=1;
    end
    return sum1+sum2
end

function algod(t,x,y)
    temp  =t[]
    for i=1 : length(x)
        push!(temp,x[i]*y[i])
    end
    sort!(temp,rev = true)     #sorting partials in reverse order
    sumd1 ::t=t(0)
    sumd2 ::t=t(0)
    j=length(temp)
    for i = 1:length(temp)
        if temp[i] < 0
            sumd1 += temp[i]  #sum from the largest to the smallest positive number for c algorithm
        end
        if temp[j]>0
            sumd2+=temp[j]
        end
        j-=1;
    end
   return sumd1+sumd2
end
x = Float32[2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y = Float32[1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

x2 = Float64[2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y2 = Float64[1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]

sum=Float32(0.0)
sum64=Float64(0)
for i=1: length(x)
    global x,y,x2,y2,sum ,sum64
    sum+=x[i]*y[i]
    sum64+=x2[i]*y2[i]
end

println("forward algorithm: 32 $(sum) and 64 $(sum64)")

sumb=Float32(0)
sumb64=Float64(0)
for i in length(x) :-1 :1
    global x,y,x2,y2,sumb,sumb64
    sumb+=x[i]*y[i]
    sumb64+=x2[i]*y2[i]
end

println("backward algorithm: 32 $(sumb) and 64 $(sumb64)")
algc32=algoc(Float32,x,y)
algc64=algoc(Float64,x2,y2)
algd32=algod(Float32,x,y)
algd64=algod(Float64,x2,y2)
println("c) $(algc32) 32 and $(algc64) 64")
println("d) $(algd32) 32 und $(algd64) 64")

