#DanielDrapala_244939
module blocksys

export  gaussian_elimination,
gaussian_elimination_with_pivots,
right_side_vector,
matrix_to_LU, solve_from_LU,
matrix_to_LU_with_pivots, solve_from_LU_with_pivots

using SparseArrays
# Funkcja obliczająca wektor prawych stron dla zadanej macierzy
# i wektora jednostkowego
# in:	A - zadana macierz, n - rozmiar macierzy
# out:	b - obliczony wektor prawych stron
function right_side_vector(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
	b = zeros(Float64, n)
	for row in 1 : n
		from_col = Int64(max(l * floor((row-1) / l) - 1, 1)) # 1 lub zaczynajac od 1 bloku gdzie znajduja sie elementy niezerowe 
		last_col = Int64(min(l + l * floor((row-1) / l), n)) # n lub konczac na macierzy A (macierz C uwzgledniona pozniej)
		for col in from_col : last_col
			b[row] += A[col, row]
		end

		if (row + l <= n)
			b[row] += A[row + l, row]   # dodanie wartosci z macierzy blokowej C 
		end
	end
	return b
end


# Funkcja rozwiązująca układ równań liniowych metodą eliminacji Gaussa
# bez wyboru elementu głównego dla macierzy o zadanej budowie
# in:	A - zadana macierz, n - rozmiar macierzy, l - wielkość bloku,
#		b - wektor prawych stron
# out:	x - rozwiązanie układu
function gaussian_elimination(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
	for pivot in 1 : n-1
		last_row = Int64(min(l + l * floor((pivot+1) / l), n))
		last_col = Int64(min(pivot + l, n)) # l to dlugosc macierzy A, k to dlugosc od pierwszej macierzy do C bez A (k= i*l+o gdzie o to indeks przekatnej macierzy C)
		for row in pivot + 1 : last_row
			if abs(A[pivot,pivot]) < eps(Float64)
				error("Współczynnik na przekątnej równy 0. Dzielenie przez zero uniemozliwia uzycie tej metody")
			end
			
			z = A[pivot, row] / A[pivot, pivot]
			A[pivot, row] = 0
			for col in pivot + 1 : last_col
				A[col, row] = A[col, row] - z * A[col, pivot]
			end
			b[row] = b[row] - z * b[pivot]
		end
	end

	x = Array{Float64}(undef,n)
	for row in n : -1 : 1
		prev_total = 0.0
		last_col = min(n, row + l)
		for col in row + 1 : last_col
			prev_total += A[col, row] * x[col]
		end
		x[row] = (b[row] - prev_total) / A[row, row]
	end
	return x
end


# Funkcja rozwiązująca układ równań liniowych metodą eliminacji Gaussa
# z częściowym wyborem elementu głównego dla macierzy o zadanej budowie
# in:	A - zadana macierz, n - rozmiar macierzy, l - wielkość bloku,
#		b - wektor prawych stron
# out:	x - rozwiązanie układu
function gaussian_elimination_with_pivots(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
	p = collect(1:n)

	for k in 1:n - 1
		last_row = Int64(min(l + l * floor((k+1) / l), n))
		last_col = Int64(min(2*l + l * floor((k+1) / l), n))
		for i in k + 1 : last_row
			max_row = k
			max = abs(A[k,p[k]])
			for x in i : last_row
				if (abs(A[k,p[x]]) > max)
					max_row = x;
					max = abs(A[k,p[x]])
				end
			end
			if (abs(max) < eps(Float64))
				error("Macierz osobliwa.")
			end
			p[k], p[max_row] = p[max_row], p[k]
			z = A[k,p[i]] / A[k,p[k]]
			A[k,p[i]] = 0
			for j in k + 1 : last_col
				A[j,p[i]] = A[j,p[i]] - z * A[j,p[k]]
			end
			b[p[i]] = b[p[i]] - z * b[p[k]]
		end
	end

	x = Array{Float64}(undef,n)
	for i in n : -1 : 1
		prev_total = 0.0
		last_col = Int64(min(2*l + l*floor((p[i]+1)/l), n))
		for j in i + 1 : last_col
			prev_total += A[j,p[i]] * x[j]
		end
		x[i] = (b[p[i]] - prev_total) / A[i, p[i]]
	end
	return x
end

# Funkcja obliczająca rozkład LU bez wyboru elementu głównego
# dla macierzy o zadanej budowie
# in:	A - zadana macierz, n - rozmiar macierzy, l - wielkość bloku
function matrix_to_LU(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
	for k in 1 : n-1
		last_row = convert(Int64, min(l + l * floor((k+1) / l), n))
		last_col = convert(Int64, min(k + l, n))
		for i in k + 1 : last_row
			if abs(A[k,k]) < eps(Float64)
				error("Współczynnik na przekątnej równy 0. Nie można zastosować metody.")
			end
			z = A[k, i] / A[k, k]
			A[k, i] = z
			for j in k + 1 : last_col
				A[j, i] = A[j, i] - z * A[j, k]
			end
		end
	end
end

# Funkcja rozwiązująca układ równań liniowych z rozkładu LU
# stworzonego bez wyboru elementu głównego
# in:	A - macierz w rozkładzie LU, n - rozmiar macierzy, l - wielkość bloku,
#		b - wektor prawych stron
# out:	x - rozwiązanie układu
function solve_from_LU(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
	z = Array{Float64}(undef,n)
	for i in 1 : n
		prev_total = 0.0
		from_col = convert(Int64, max(l * floor((i-1) / l) - 1, 1))
		for j in from_col : i-1
			prev_total += A[j, i] * z[j]
		end
		z[i] = b[i] - prev_total
	end

	x = Array{Float64}(undef,n)
	for i in n : -1 : 1
		prev_total = 0.0
		last_col = min(n, i + l)
		for j in i + 1 : last_col
			prev_total += A[j, i] * x[j]
		end
		x[i] = (z[i] - prev_total) / A[i, i]
	end
	return x
end

# Funkcja obliczająca rozkład LU z częściowym wyborem elementu głównego
# dla macierzy o zadanej budowie
# in:	A - zadana macierz, n - rozmiar macierzy, l - wielkość bloku
# out:	p - wektor permutacji wierszy
function matrix_to_LU_with_pivots(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
	p = collect(1:n)

	for k in 1:n - 1
		last_row = Int64(min(l + l * floor((k+1) / l), n))
		last_col = Int64(min(2*l + l*floor((k+1)/l), n))
		for i in k + 1 : last_row
			max_row = k
			max = abs(A[k,p[k]])
			for x in i : last_row
				if (abs(A[k,p[x]]) > max)
					max_row = x;
					max = abs(A[k,p[x]])
				end
			end
			if (abs(max) < eps(Float64))
				error("Macierz osobliwa.")
			end
			p[k], p[max_row] = p[max_row], p[k]
			z = A[k,p[i]] / A[k,p[k]]
			A[k,p[i]] = z
			for j in k + 1 : last_col
				A[j,p[i]] = A[j,p[i]] - z * A[j,p[k]]
			end
		end
	end
	return p
end

# Funkcja rozwiązująca układ równań liniowych z rozkładu LU
# stworzonego z częściowym wyborem elementu głównego
# in:	A - macierz w rozkładzie LU, n - rozmiar macierzy, l - wielkość bloku,
#		b - wektor prawych stron, p - wektor permutacji wierszy
# out:	x - rozwiązanie układu
function solve_from_LU_with_pivots(A::SparseMatrixCSC{Float64, Int64}, n::Int64,
							l::Int64, b::Vector{Float64}, p::Vector{Int64})
	z = Array{Float64}(undef,n)
	for i in 1 : n
		prev_total = 0.0
		from_col = Int64(max(l * floor((p[i]-1) / l) - 1, 1))
		for j in from_col : i-1
			prev_total += A[j, p[i]] * z[j]
		end
		z[i] = b[p[i]] - prev_total
	end

	x = Array{Float64}(undef,n)
	for i in n : -1 : 1
		prev_total = 0.0
		last_col = Int64(min(2*l + l*floor((p[i]+1)/l), n))
		for j in i + 1 : last_col
			prev_total += A[j, p[i]] * x[j]
		end
		x[i] = (z[i] - prev_total) / A[i, p[i]]
	end
	return x
end
end