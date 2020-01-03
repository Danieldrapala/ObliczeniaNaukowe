#DanielDrapala_244939


module fileinout
export load_matrix_ffile, load_vector_ffile, save_values

using SparseArrays,LinearAlgebra
# wczytywanie macierzy z pliku 
# pierwsze wartosci n l potem i j A[i,j]
# in:	f - ścieżka do pliku z macierzą
# out:	A - macierz rzadka wczytana z pliku, n - rozmiar macierzy,
#		l - wielkość bloku
function load_matrix_ffile(f::String)
	open(f) do file
		ln = split(readline(file))
		n = parse(Int64, ln[1])
		l = parse(Int64, ln[2])
		el_num = n*l + 3*(n-l) 
		J = Array{Int64}(undef,el_num)
		I = Array{Int64}(undef,el_num)
		V = Array{Float64}(undef,el_num)
		it = 1
		while !eof(file)
			ln = split(readline(file))
			J[it] = parse(Int64, ln[1])
			I[it] = parse(Int64, ln[2])
			V[it] = parse(Float64, ln[3])
			it += 1
		end
		A = sparse(I, J, V)   #I- kolumny J-wiersze V-wartosci
		return (A, n, l)
	end
end


# wczytywanie wektoru prawych stron z pliku
# in:	f - ścieżka do pliku z wektorem
# out:	b - wczytany wektor
function load_vector_ffile(f::String)
	open(f) do file
		n = parse(Int64, readline(file))
		b = Array{Float64}(undef,n)
		it = 0
		while !eof(file)
			it += 1
			b[it] = parse(Float64, readline(file))
		end
		return b
	end
end


# Funkcja zapisująca rozwiązanie układu równań do pliku
# in:	f - ścieżka do pliku zapisowego, x - wektor z rozwiązaniem,
#		n - rozmiar macierzy, err - wybór zapisu błędu względnego
function save_values(f::String, x::Array{Float64}, n::Int64, err::Bool)
	open(f, "w") do file
		if (err)
			relative_error = norm(ones(n) - x) / norm(x)
			println(file, relative_error)
		end
		for i in 1 : n
			println(file, x[i])
		end
	end
end
end