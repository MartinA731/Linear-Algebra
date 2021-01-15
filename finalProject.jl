#Matrix that you want to find (n x n)
A = [ 0 1 ; 1 0]

#initial guess vector, must have n components and be nonzero
#This vector must have at least one component in the direction
#of the largest eigenvector you are looking for
#recommended that you test with several guesses to guarantee 
#there are no errors
v0 = [ 2 ; 2 ]

#The highest possible error value for powerMethod function below
bound = 10^-35

keepIterating = true

#power method function that returns eigenvectors that correspond
#to the highest eigenvalue
function powerMethod(A, x, keepIterating, tolerable_error)
    #placeholders
    lambda_old = 1
    eigenvalue = 1
    eigenvector = 1
    #Continue iteration as long as vectors don't converge
    while(keepIterating)
        #case for the 0 eigenvalue
        sumElements = 0
        test = A*x
        for i in 1 : size(test)[1]
            sumElements += test[i]^2
        end
        if(sumElements == 0)
            return [0, x]
        end
        #continue iterating through new possible eigenvalues and
        #eigenvectors until the eigenvalues and vectors converge
        x = A*x
        lambda_new = maximum(abs.(x)) 
        if(lambda_new != 0)
            x = x/lambda_new
        end
        error = abs(lambda_new - lambda_old)
        keepIterating = error > tolerable_error 
        lambda_old = lambda_new
        if(keepIterating == false)
            eigenvalue = lambda_new
            eigenvector = x
        end
    end
    return eigenvector
end

max_eigenvector = powerMethod(A,v0,true,bound)

#Use the rayleigh quotient to find eigenvalue that corresponds
#to eigenvalue from power method
function rayleighQuotient(A, eigenvec_1)
    return (A*eigenvec_1)'*eigenvec_1 / (eigenvec_1'*eigenvec_1)
end

max_lambda = rayleighQuotient(A,max_eigenvector)

#Make the first eigenvector that we calculated a unit 
#eigenvector for correct results for altered matrix B
sumElements = 0
length = size(max_eigenvector)[1]
for i in 1 : length
    sumElements += max_eigenvector[i]^2
end
scaleUnit = 1 / sumElements^(1/2)
unit_eigenvector = scaleUnit*max_eigenvector

#create matrix B with the same same eigenvectors and
#eigenvalues except the largest one has been replaced by 0. 
#We can use this to calculate second highest eigenvalue,eigenvector
B = A - max_lambda*unit_eigenvector*unit_eigenvector'

#we have to change our guess vector because (2, 2) creates 
#errors for matrix B
v0 = [7 ; 20]

second_eigenvector = powerMethod(B,v0,true,bound)

second_eigenvalue = rayleighQuotient(A,second_eigenvector)
