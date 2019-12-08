include("metody.jl")
using  Main.metody

f(x) = x 
g(x) = x^2 

rysujNnfx(f, -10.0, 10.0, 3)
rysujNnfx(f, -10.0, 10.0, 20)

rysujNnfx(g, -10.0, 10.0, 3)
rysujNnfx(g, -10.0, 10.0, 20)
