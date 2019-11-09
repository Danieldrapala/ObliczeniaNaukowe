#Daniel Drapala 244939
using Polynomials

P= Float64[1, -210.0-2.0^(-23.0), 20615.0,-1256850.0,
      53327946.0,-1672280820.0, 40171771630.0, -756111184500.0,          
      11310276995381.0, -135585182899530.0,
      1307535010540395.0,     -10142299865511450.0,
      63030812099294896.0,     -311333643161390640.0,
      1206647803780373360.0,     -3599979517947607200.0,
      8037811822645051776.0,      -12870931245150988800.0,
      13803759753640704000.0,      -8752948036761600000.0,
      2432902008176640000.0]

P=P[end:-1:1]

p = Float64[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
      14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]

P_Poly = Poly(P)
p_poly = poly(p)
calculated_roots = roots(P_Poly)

for i = 1:length(calculated_roots)
      println("$(calculated_roots[i] )")
end
for i = 1:length(calculated_roots)
    println("$(abs(polyval(P_Poly,calculated_roots[i])))")
end
for i = 1:length(calculated_roots)
      println("$(abs(polyval(p_poly,calculated_roots[i])))")
end
for i = 1:length(calculated_roots)
      println("$(abs(calculated_roots[i] - i))")
end