function [ x, err ] = gradient( x, A, b, iterations )
%GRADIENT Does gradient descent relaxation.
%   x - Initial guess for solution to Ax=b.
%   A - A of Ax=b
%   b - b of Ax=b
% iterations - How many relaxation iterations to do

err = zeros([iterations,1]);

for i=1:iterations
    r = b-A*x;
    lambda = r'*r/(r'*A*r);
    x = x + lambda*r;
    err(i) = norm(A*x-b);
    fprintf(1,'Iteration %d: %.3f\n', i, err(i));
end

end

