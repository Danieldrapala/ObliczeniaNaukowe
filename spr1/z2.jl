#author: Daniel Drapala 244939
# Kahan's machine epsilon

function kahansmethod(t)
    keps::t = t(3.0)*(t(4.0)/t(3.0) - t(1.0)) - t(1.0)
    return keps
end

for t in [Float16, Float32, Float64]
    println("Calculated Kahan's eps for $(t): $(kahanseps(t))")
    println("      Function macheps for $(t): $(eps(t))\n")
end