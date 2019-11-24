include("metody.jl")
using  Main.metody
using Test
delta = 0.5 * 10.0^-5.0             
epsilon = 0.5 * 10.0^-5.0 

f(x) = 0.5*(x)^2.0 -8
pf(x)= 2*x
g(x) = x^2 - 16
pg(x) = 2*x
j(x) = x^4 - 16
pj(x) = 4*(x^3)
(r,v,it,err)=mbisekcji(f,Float64(1),Float64(10),delta,epsilon)
(r1,v1,it1,err1)=mstycznych(f,pf,Float64(1),delta,epsilon,40)
(r2,v2,it2,err2)=msiecznych(f,Float64(1),Float64(3),delta,epsilon,40)

println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
println("r: $(r1)\t v: $(v1)\t it: $(it1)\t err: $(err1)")
println("r: $(r2)\t v: $(v2)\t it: $(it2)\t err: $(err2)")

@testset "Metods polynomial test" begin
@test ((r <= 4.0 + epsilon) & (r >= 4.0-epsilon))
@test ((r1 <= 4.0 + epsilon) & (r1 >= 4.0-epsilon))
@test ((r2 <= 4.0 + epsilon) & (r2 >= 4.0-epsilon))

end
(r,v,it,err)=mbisekcji(g,Float64(1),Float64(10),delta,epsilon)
(r1,v1,it1,err1)=mstycznych(g,pg,Float64(1),delta,epsilon,40)
(r2,v2,it2,err2)=msiecznych(g,Float64(1),Float64(2),delta,epsilon,40)

println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
println("r: $(r1)\t v: $(v1)\t it: $(it1)\t err: $(err1)")
println("r: $(r2)\t v: $(v2)\t it: $(it2)\t err: $(err2)")

@testset "Metods polynomial test" begin
@test ((r <= 4.0 + epsilon) & (r >= 4.0-epsilon))
@test ((r1 <= 4.0 + epsilon) & (r1 >= 4.0-epsilon))
@test ((r2 <= 4.0 + epsilon) & (r2 >= 4.0-epsilon))

end
(r,v,it,err)=mbisekcji(j,Float64(1),Float64(4),delta,epsilon)
(r1,v1,it1,err1)=mstycznych(j,pj,Float64(1),delta,epsilon,40)
(r2,v2,it2,err2)=msiecznych(j,Float64(1),Float64(1.5),delta,epsilon,40)

println("r: $(r)\t v: $(v)\t it: $(it)\t err: $(err)")
println("r: $(r1)\t v: $(v1)\t it: $(it1)\t err: $(err1)")
println("r: $(r2)\t v: $(v2)\t it: $(it2)\t err: $(err2)")

@testset "Metods polynomial test" begin
    
    @test ((r <= 2.0 + epsilon) & (r >= 2.0-epsilon))
    @test ((r1 <= 2.0 + epsilon) & (r1 >= 2.0-epsilon))
    @test ((r2 <= 2.0 + epsilon) & (r2 >= 2.0-epsilon))

end