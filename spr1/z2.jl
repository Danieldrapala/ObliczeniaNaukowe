#author: Daniel Drapala 244939
# Kahan's machine epsilon

using Printf
function kahansmethod(t)
    keps::t = t(3.0)*(t(4.0)/t(3.0) - t(1.0)) - t(1.0)
    return keps
end

for t in [Float16, Float32, Float64]
    @printf("Calculated Kahan's eps for %s : %.20e\n",t,kahansmethod(t))
    @printf("Function macheps for %s : %.20e\n",t,eps(t))
end