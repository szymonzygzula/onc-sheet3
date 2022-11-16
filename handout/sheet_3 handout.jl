# using IterativeSolvers

include("sheet_3_utility.jl")

###############
# Exercise 2a #
###############

function smoothing(A, b, x, maxiter = 3)
    res = zeros(size(b)[1])
    return res
end

###############
# Exercise 2b #
###############

function residual(A, b, x)
    return zeros(size(b)[1])
end

###############
# Exercise 2c #
###############
function computeInterpolationMatrix(m, n)
    res = zeros(m,n)
    return res
end

function computeRestrictionMatrix(m, n)
    res = zeros(m,n)
    return res
end


function restriction(A, r)

    # restricted matrix
    A_r = zeros(size(A)[1], size(A)[2])

    # restricted residual
    b_r = zeros(size(r)[1])

    return A_r, b_r
end

###############
# Exercise 2d #
###############

function prolongation(x_e)
    return x_e
end

function V_Cycle(A, b, step=1)

    # Smoothing
    x_solution = smoothing(A, b, x);
    # x_solution = gauss_seidel(A,b,maxiter=3)

    # Compute Residual Errors
    r = residual(A, b, x_solution);

    # Restriction
    A_r, b_r = restriction(A, r);

    # Stop recursion at smallest grid size, otherwise continue recursion
    if true
        x_e = zeros(size(b_r)[1])
        x_e = smoothing(A_r, b_r, x_e);
        # x_e = gauss_seidel(A_r, b_r, maxiter=3);
    else
        x_e = V_Cycle(A_r, b_r, step+1);
    end

    # Compute the interpolated solution from the residual
    x_e_p = prolongation(x_e);
    
    # Correction
    x_solution += x_e_p;

    # Smoothing
    x_solution = smoothing(A, b, x_solution);
    # x_solution = gauss_seidel!(x_solution,A,b,maxiter=3)

    return x_solution
end

###############
# Exercise 3a #
###############

function F_Cycle(A, b)

end

###############
# Exercise 3b #
###############

function W_Cycle(A, b)

end

###############
# Exercise 2a #
###############

function jacobi_solver(A::Matrix, b::Vector, x_start::Vector, iter::Int64)
    x_prev = x_start 
    x_next = deepcopy(x_prev)
    
    error = Float64[]
    
    D = Diagonal(A)
    L = LowerTriangular(A) - D
    U = UpperTriangular(A) - D
    D_inv = inv(D)

    for k = 1:iter
        x_next = D_inv * (b - (L + U) * x_prev)

        push!(error, norm(x_next - x_prev))
        
        x_prev = deepcopy(x_next)
    end

    x_next, error
end

m = 10

A = create1DPoissonMatrix(m)
b = createCoefficientVector(m)

x, r = iterative_solver(A, b)

println(Matrix{Float64}(I, m, m) * length(x))

r_H = restriction(r)

x = zeros(m)
multigrid_result = V_Cycle(A,b,x)

# Initializing linear system data for exercise 2
A = [6.0 1.0 2.0; 3.0 5.0 1.0; 1.0 1.0 4.0]
x = [1.0, 0.0, 1.0]
b = [1.0, 2.0, 3.0]

# Calculating result matrix x and error (5 iterations) for exercises 2a and 2b
iter_num = 5
x, e = jacobi_solver(A, b, x, iter_num)

# Priting results for excerice 2
# println("2a) Solution x: ", x)
# println("2b) Calculated error (only last iteration): ", e[end])
# println("2b) Calculated error (all iterations): ", e)   