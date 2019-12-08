#author: Daniel Drapala

module metody
export ilorazyRoznicowe, warNewton, naturalna, rysujNnfx
using Plots
plotly()

#ZADANIE 1
# x- wektor wezlow
function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    
        n = length(f)             # długość wektorów
        fx = Vector{Float64}(undef,n)   # deklaracja wektora z ilorazami różnicowymi
    
        for i = 1 : n             # przekopiowanie wartości funkcji interpolowanej
            fx[i] = f[i]
        end
    
        for i = 2 : n              # obliczanie ilorazów różnicowych
            for j = n : -1 : i
                fx[j] = (fx[j] - fx[j - 1]) / (x[j] - x[j - i + 1])
            end
        end
    
        return fx
    end

# ZADANIE 2

# Uogólnionego algorytm Hornera
# x - wektor zawierający węzły, 
# fx - wektor zawierający ilorazy różnicowe,
# arg - punkt w którym należy obliczyć wartość wielomianu

# out: val - wartość wielomiau w punkcie arg

function warNewton(x :: Vector{Float64}, fx :: Vector{Float64}, arg :: Float64)
    val = Float64(1.0)            # deklaracja outputu
    n = length(fx)                # długość wektorów
    val = fx[n]
	for i = n-1 : -1 : 1
		val = fx[i]+(arg-x[i])*val
    end
    
    return val
end


#ZADANIE 3

# Funkcja obliczająca postac normalna wspolczynnikow posiadajac wielomain interpolacyjny Newtona
# x  - wektor zawierający węzły,
# fx - wektor zawierający ilorazy różnicowe
# a  - wektor zawierający współczynniki w postaci normalnej
function naturalna(x :: Vector{Float64}, fx :: Vector{Float64})
    n = length(fx)                  
    a = Vector{Float64}(undef,n)          # współczynniki w postaci normalnej
    a[n] = fx[n]                          # z twierdzenia a_n = c_n
    for k = n-1 : -1 : 1            
        a[k] = fx[k]-a[k+1]*x[k]    
        for i = k+1 : n-1           
            a[i] = a[i]-a[i+1]*x[k]
        end
    end
    return a
end


# ZADANIE 4

# F rysuje funkcje i jej interpertacje interpolacyjna
# f -dana funkcja; 
# a, b - końce przedziału; 
# n - stopień wielomianu

function rysujNnfx(f, a :: Float64, b :: Float64, n :: Int)
    x = Vector{Float64}(undef,n+1)        # wektory potrzebne do stworzenia wielomianu 
    y = Vector{Float64}(undef,n+1)        
    fx = Vector{Float64}(undef,n+1)       

    acc = 30                       # mnożnik dla dokładiejszego rysowania wykresów

    ploty = Vector{Float64}(undef,acc*(n + 1))       # wektory do wykresow powiekszone acc razy
    plotx = Vector{Float64}(undef,acc*(n + 1))      
    plotinter = Vector{Float64}(undef,acc*(n + 1))    

    kh = 0.0                                   
    max_n = n + 1                              
    h = (b-a)/n                                

    for i = 1: max_n
        x[i] = a + kh
        y[i] = f(x[i])
        kh += h
    end
    fx = ilorazyRoznicowe(x, y);
    max_n *= acc
    kh = 0.0
    h = (b - a)/(max_n-1)

    for i = 1: max_n
        plotx[i] = a + kh
        plotinter[i] = warNewton(x, fx, plotx[i])
        ploty[i] = f(plotx[i])
        kh += h
    end
    titles=string("n=",n)
    plt = plot(plotx, [ploty, plotinter],title=titles, color= [:green :red],label=["funkcja" "interpolacja"], linewidth=2.0)
    display(plt)
end
end