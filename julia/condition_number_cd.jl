# Import required modules
using LinearAlgebra, SparseArrays, DelimitedFiles

path = "../output_files/paper/example2/data/"

n = 6
condition_numbers = zeros(n)

for i=1:n
    println(string(i) * "/" * string(n))
    A = readdlm(path * "mat_" * string(i) * ".dat")
    A = sparse(A[:,1], A[:,2], A[:,3])

    #condition_numbers[i] = cond(A, 1)           # condest norm
    condition_numbers[i] = cond(Matrix(A), 2)   # spectral norm
    println(string(condition_numbers[i]) * "\n")
end

# Display the result
println("condition_number = ", condition_numbers)