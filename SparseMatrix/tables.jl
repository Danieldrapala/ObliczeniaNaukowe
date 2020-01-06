include("blocksys.jl")
include("matrixgen.jl")
include("fileinout.jl")

using Main.blocksys
using Main.matrixgen
using Main.fileinout
using LinearAlgebra
using SparseArrays
block_size = 4
gen_sizes = [52, 500, 5000, 25000,50000]
# 
function compare_all()
    for size in gen_sizes
        blockmat(size, block_size, 25.0, "matgenerated/$(size)/A.txt")
        A, n, l = load_matrix_ffile("matgenerated/$(size)/A.txt")
        b = right_side_vector(A, size, block_size)
        x = ones(Float64, n)
        A1, b1 = deepcopy(A), deepcopy(b)
        Ap, bp = deepcopy(A), deepcopy(b)
        Al, bl = deepcopy(A), deepcopy(b)
        Alu, blu = deepcopy(A), deepcopy(b)
       

        # if (size<25001) 
        #     Atypical, btypical = Array(A), deepcopy(b)
        #     result = @timed \(Atypical, btypical)
        # end
        gaussian_elimination(A1,n,l,b1,false)
        dropzeros(A1)
        
        ge = @timed solving_eq_after_GE(A1, n, l, b1)
        
        # gewp=  @timed gaussian_elimination_with_pivots(Ap,n,l,bp,false)
        gaussian_elimination_with_pivots(Ap,n,l,bp,false)
        dropzeros(Ap)
        
        gewp = @timed solving_eq_after_GEWP(Ap, n, l, bp)
        
        gaussian_elimination_with_pivots(Al,n,l,bl,true)
        dropzeros(Al) 
        # lu= @timed gaussian_elimination(Al,n,l,bl,true)
        lu = @timed solving_eq_after_LU(Al, n, l, bl)
       

        # luwp=  @timed gaussian_elimination_with_pivots(Alu,n,l,blu,true)
        gaussian_elimination_with_pivots(Alu,n,l,blu,true)
        dropzeros(Alu)

        luwp = @timed solving_eq_after_LUWP(Alu, n, l, blu) 
        # println("$n :"*(size<250001 ?  "$(norm(x - result[1]) / norm(result[1])):" : "")* " $(norm(x - ge[1]) / norm(ge[1])) : $(norm(x - gewp[1]) / norm(gewp[1])) : $(norm(x - lu[1]) / norm(lu[1])) : $(norm(x - luwp[1]) / norm(luwp[1]))\\\\")
        # * (size<250001 ? "$(round(result[2],digits=6)) : $(round(result[3]/2^20, digits=4)):" : "")
        println("$size :"*"$(round(ge[2], digits=6)) : $(round(ge[3]/2^20, digits=4)) : $(round(gewp[2], digits=6)) : $(round(gewp[3]/2^20, digits=4)) : $(round(lu[2], digits=6)) : $(round(lu[3]/2^20, digits=4)) : $(round(luwp[2], digits=6)) : $(round(luwp[3]/2^20, digits=4)) \\\\")
     end

end


function gaussian_eliminationhehe(A::Array{Float64,2}, n::Int64, l::Int64, b::Vector{Float64},LU::Bool)
    for pivot in 1 : n-1
        			
        for row in pivot + 1 : n						
            if abs(A[pivot,pivot]) < eps(Float64)
                error("Współczynnik na przekątnej równy 0. Dzielenie przez zero uniemozliwia uzycie tej metody")
            end
            
            z = A[row,pivot] / A[pivot, pivot]
            A[row,pivot] = LU ? z : 0

            for col in pivot + 1 : n
                A[row,col] = A[row,col] - z * A[pivot,  col]
            end
            if !LU
                b[row] = b[row] - z * b[pivot]
            end
        end
    end

    x=zeros(n)
    
    
    for row in n : -1 : 1
        prev_total = 0.0
        
        for col in row + 1 : n
            prev_total += A[row,col ] * x[col]
        end
        x[row] = (b[row] - prev_total) / A[row, row]
    end
    return x
end

compare_all()