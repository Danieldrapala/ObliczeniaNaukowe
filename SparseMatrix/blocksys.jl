#DanielDrapala_244939
module blocksys

export  
right_side_vector,
solving_eq_after_GE, solving_eq_after_GEWP,
solving_eq_after_LU, solving_eq_after_LUWP,gaussian_elimination_with_pivots, solving_eq_after_GE_with_no_matrix,gaussian_elimination

using SparseArrays
# Funkcja obliczająca wektor prawych stron dla zadanej macierzy
# i wektora jednostkowego
# wejscie:	A - zadana macierz, 
#			n - rozmiar macierzy
# wyjscie:	b - obliczony wektor prawych stron
	function right_side_vector(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64)
		b = zeros(Float64, n)
		for row in 1 : n
			from_col = Int(max(l * (floor((row-1) / l) - 1), 1)) 	# 1 lub zaczynajac od 1 bloku gdzie znajduja sie elementy niezerowe 
			last_col = Int(min(l + l * floor((row-1) / l), n)) 	# n lub konczac na macierzy A (macierz C uwzgledniona pozniej)
			for col in from_col : last_col
				b[row] += A[row,col]
			end

			if (row + l <= n)
				b[row] += A[row , row + l]  				 # dodanie wartosci z macierzy blokowej C 
			end
		end
		return b
	end


	# Funkcja stosujaca metode eliminacji Gaussa
	# bez wyboru elementu głównego dla macierzy o zadanej budowie
	# z flaga LU  zapisuje rowniez macierz w rozkladzie LU 
	# wejscie:	
	#			A - zadana macierz, 
	#			n - rozmiar macierzy, 
	#			l - wielkość bloku,
	#			b - wektor prawych stron
	#			LU- flaga zapisu macierzy w rozkladzie LU
	function gaussian_elimination(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64},LU::Bool)
		for pivot in 1 : n-1
			last_row = Int(min(l + l * floor((pivot+1) / l), n))# l to dlugosc macierzy A, k to dlugosc od pierwszej macierzy do C bez A 
			last_col = min(pivot + l, n)				
			for row in pivot + 1 : last_row						
				if abs(A[pivot,pivot]) < eps(Float64)
					error("Współczynnik na przekątnej równy 0. Dzielenie przez zero uniemozliwia uzycie tej metody")
				end
				
				z = A[row,pivot] / A[pivot, pivot]
				A[row,pivot] = LU ? z : 0

				for col in pivot + 1 : last_col
					A[row,col] = A[row,col] - z * A[pivot,  col]
				end
				if !LU
					b[row] = b[row] - z * b[pivot]
				end
			end
		end
	end
	#Funkcja rozwiazujaca uklad rownan w macierzy uprzednio uporzadkowanej eliminacja Gaussa
	#solving_eq_after_gaussian_elimination
	
	function solving_eq_after_GE(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		x=zeros(n)
		gaussian_elimination(A,n,l,b,false)
		
		for row in n : -1 : 1
			prev_total = 0.0
			last_col = Int(min(n, row + l))
			for col in row + 1 : last_col
				prev_total += A[row,col ] * x[col]
			end
			x[row] = (b[row] - prev_total) / A[row, row]
		end
		return x
	end
	function solving_eq_after_GE_with_no_matrix(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		x=zeros(n)
		
		for row in n : -1 : 1
			prev_total = 0.0
			last_col = Int(min(n, row + l))
			for col in row + 1 : last_col
				prev_total += A[row,col ] * x[col]
			end
			x[row] = (b[row] - prev_total) / A[row, row]
		end
		return x
	end
	# Funkcja przeksztalcajaca macierz metodą eliminacji Gaussa
	# z częściowym wyborem elementu głównego dla macierzy o zadanej budowie
	# wejscie:	
	#		A - zadana macierz, 
	#		n - rozmiar macierzy, 
	#		l - wielkość bloku,
	#		b - wektor prawych stron
	function gaussian_elimination_with_pivots(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64},LU::Bool)
		p = collect(1:n)

		for pivot in 1:n - 1
			last_row = Int(min(l + l * floor((pivot+1) / l), n))
			last_col = Int(min(2*l + l * floor((pivot+1) / l), n))
			for row in pivot + 1 : last_row
					max_row = pivot
					max = abs(A[p[pivot],pivot])

				for i in row : last_row
					if (abs(A[p[i],pivot]) > max)
						max_row,max = i,abs(A[p[i],pivot]);
					end
				end

				if (abs(max) < eps(Float64))
					error("Macierz osobliwa.")
				end
				p[pivot], p[max_row] = p[max_row], p[pivot]
				
				z = A[p[row],pivot] / A[p[pivot],pivot]

				A[p[row],pivot] = LU ? z : 0		

				for col in pivot + 1 : last_col
					A[p[row],col] = A[p[row],col] - z * A[p[pivot],col]
				end
				if !LU
					b[p[row]] = b[p[row]] - z * b[p[pivot]]
				end
			end
		end
			return p
		
	end
	#solving_eq_after_gaussian_elimination_with_pivots
	function solving_eq_after_GEWP(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		# Perform gaussian elimination
		p=gaussian_elimination_with_pivots(A, n, l, b, false)

		x = zeros(n)
		for row in n : -1 : 1
			prev_total = 0
			last_col = Int(min(2*l + l*floor((p[row]+1)/l), n))
			for col in row + 1 : last_col
				prev_total += A[p[row],col] * x[col]
			end
			x[row] = (b[p[row]] - prev_total) / A[p[row],row]
		end
		return x
	end
	function solving_eq_after_GEWP_nogewp(p::Vector{Int64},A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		x = zeros(n)
		for row in n : -1 : 1
			prev_total = 0
			last_col = Int(min(2*l + l*floor((p[row]+1)/l), n))
			for col in row + 1 : last_col
				prev_total += A[p[row],col] * x[col]
			end
			x[row] = (b[p[row]] - prev_total) / A[p[row],row]
		end
		return x
	end

	function LU(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
			return gaussian_elimination(A,n,l,b,true)
	end

	function LUWP(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		return gaussian_elimination_with_pivots(A,n,l,b,true)
	end

# Funkcja rozwiązująca układ równań liniowych z rozkładu LU
# stworzonego bez wyboru elementu głównego
# in:	A - macierz w rozkładzie LU, n - rozmiar macierzy, l - wielkość bloku,
#		b - wektor prawych stron
# out:	x - rozwiązanie układu
	function solving_eq_after_LU(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		LU(A,n,l,b)
		z = zeros(n)
		for row in 1 : n
			prev_total = 0.0
			from_col = Int(max(l * floor((row-1) / l) - 1, 1))
			for col in from_col : row-1
				prev_total += A[row,col] * z[col]
			end
			z[row] = b[row] - prev_total
		end

		x = zeros(n)
		for row in n : -1 : 1
			prev_total = 0.0
			last_col = min(n, row + l)
			for col in row + 1 : last_col
				prev_total += A[row, col] * x[col]
			end
			x[row] = (z[row] - prev_total) / A[row, row]
		end
		return x

	end

# Funkcja rozwiązująca układ równań liniowych z rozkładu LU
# stworzonego z częściowym wyborem elementu głównego
# wejscie:	A - macierz w rozkładzie LU, 
#			n - rozmiar macierzy, l - wielkość bloku,
#			b - wektor prawych stron, p - wektor permutacji wierszy
# wyjscie:	x - rozwiązanie układu
	function solving_eq_after_LUWP(A::SparseMatrixCSC{Float64, Int64}, n::Int64, l::Int64, b::Vector{Float64})
		p=LUWP(A,n,l,b)
		z = zeros(n)
		for row in 1 : n
			prev_total = 0.0
			from_col = Int(max(l * floor((p[row]-1) / l) - 1, 1))
			for col in from_col : row-1
				prev_total += A[ p[row],col] * z[col]
			end
			z[row] = b[p[row]] - prev_total
		end

		x=zeros(n)
		for row in n : -1 : 1
			prev_total = 0.0
			last_col = Int(min(2*l + l*floor((p[row]+1)/l), n))
			for col in row + 1 : last_col
				prev_total += A[p[row],col] * x[col]
			end
			x[row] = (z[row] - prev_total) / A[p[row], row]
		end
		return x
	end

end