using LinearAlgebra

###############
# Exercise 1a 
###############
function createCoefficientVector(m)
    res = zeros(m)

    for i = 1:m
        res[i] = rand(Float64) * 5  
    end

    return res
end

###############
# Exercise 1b #
###############

function create1DPoissonMatrix(m)
    res = zeros(m, m)

    for i = 1:m
        for j = 1:m
            if i == j
                res[i, j] = 2
            elseif abs(i - j) == 1
                res[i, j] = -1
            end
        end
    end

    res
end

###############
# Exercise 1c #
###############

function iterative_solver(A::Matrix, b::Vector, maxiter = 100, epsilon = 1e-8)
	x_prev = zeros(length(b))  
    x_next = deepcopy(x_prev)

    for k = 1:maxiter
        for i = 1:length(b)
            sum = 0
            for j = 1:length(b)
                if i != j
                    sum += A[i, j] * x_prev[j]
                end
            end
    
            x_next[i] = 1 / A[i, i] * (b[i] - sum)
        end
        
		if norm(b - A * x_next) <= epsilon
			break
        end
        
        x_prev = deepcopy(x_next)
        
    end

    x_next, A * b - x_next
end

###############
# Exercise 1d #
###############

function restriction(r)
    r_H = zeros(length(r))
    

    r_H
end