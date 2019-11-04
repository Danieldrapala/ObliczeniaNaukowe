#author Daniel Drapala 244939

using Printf

#m- machine epsilon for type t (float 16,32 or 64)
function macheps(t)
    m::t = t(1.0)
    while t(1.0) + m/t(2.0) > t(1.0)
        m /= t(2.0)
    end
    return m
end

#eta- the smallest number eta >0 for type t (float 16,32 or 64)
function eta(t)
    eta::t = t(1.0)
    while eta/t(2.0) > t(0.0)
        eta /= t(2.0)
    end
    return eta
end

#max- the maximum number for each  type t(float 16,32,64)
function flomax(t)
    max::t = t(1.0)
    while !isinf(max*t(2.0))
        max *= t(2.0)
    end
    max*=(t(2.0)-t(macheps(t)))
    return max
end


for t in [Float16, Float32, Float64]
    @printf("Calculated macheps for %s: %.15e\n",t, macheps(t))
    @printf("Function   macheps for %s: %.15e\n\n",t, eps(t))

    @printf("Calculated eta for %s: %.15e\n",t, eta(t))
    @printf("Function   eta for %s: %.15e\n\n",t, nextfloat(t(0.0)))

    @printf("Calculated max for %s: %.15e\n",t, flomax(t))
    @printf("Function   max for %s: %.15e\n\n",t, floatmax(t))
end