include("metody.jl")
using Main.metody

# Funkcja e^(1-x)-1, x - argument
f(x) = MathConstants.e^(1.0 - x) - 1.0
# Funkcja xe^-x, x - argument
g(x) = x * MathConstants.e^(-x)

# Pochodna funkcji f
pf(x) = -MathConstants.e^(1.0 - x)
# Pochodna funkcji g
pg(x) = MathConstants.e^-x - x *  MathConstants.e^-x

delta = 10.0^-5.0               # dokładność obliczeń
epsilon = 10.0^-5.0             # dokładność obliczeń
maxit = 500                      # maksymalna dopuszczalna liczba iteracji

println("y = e^(1-x)-1: ")

# # println("Metoda bisekcji: ")
# # (r,v,it,err) = mbisekcji(f, 0.0, 1.5, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(f, 0.5, 3.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(f, 0.0, 4.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(f, 0.0, 50.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")

# # println("Metoda stycznych: ")
# # (r,v,it,err) = mstycznych(f, pf, -1.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(f, pf, 0.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(f, pf, 1.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(f, pf, 5.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(f, pf, 7.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# (r,v,it,err) = mstycznych(f, pf, 8.0, delta, epsilon, 1000000000)
# println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(f, pf, 20.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")


println("Metoda siecznych: ")
(r,v,it,err) = msiecznych(f, -1.0, 2.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(f, 0.5, 2.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(f, -3.0, 6.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(f, -2.0, 10.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(f, 10.0, 50.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")

# println("y = xe^-x: ")

# # println("Metoda bisekcji: ")
# # (r,v,it,err) = mbisekcji(g, 0.0, 1.5, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(g, 0.5, 3.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(g, 0.0, 4.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mbisekcji(g, 0.0, 50.0, delta, epsilon)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")

# # println("Metoda stycznych: ")
# # (r,v,it,err) = mstycznych(g, pg, -1.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 0.5, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 1.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 5.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 7.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 8.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
# # (r,v,it,err) = mstycznych(g, pg, 20.0, delta, epsilon, maxit)
# # println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
println("Metoda siecznych: ")
(r,v,it,err) = msiecznych(g, -1.0, 0.5, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(g, -0.25, 1.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(g, -3.0, 1.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(g, 10.0, 20.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
(r,v,it,err) = msiecznych(g, -2.0, 50.0, delta, epsilon, maxit)
println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")

