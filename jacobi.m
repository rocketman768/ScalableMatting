function [ x, err ] = jacobi( x, A, b, iterations, omega )
%JACOBI Does (damped) Jacobi relaxation.
%   x - Initial guess for solution to Ax=b.
%   A - A of Ax=b
%   b - b of Ax=b
% iterations - How many Jacobi relaxation iterations to do
% omega - Damping factor [0,1]. 1 means no damping, and 0 means infinite.

err = zeros([iterations,1]);

% This is traditional Jacobi relaxation.
d = diag(A,0);
bias = b./d;

% This is a more stable relaxation.
%d = max(diag(A,0));
%bias = b./d;

for i=1:iterations
    x = (x - (A*x)./d + bias)*omega + x*(1-omega);
    %imshow(reshape(x,[552,800]),[0,1]); pause(0.2);
    err(i) = norm(A*x-b);
    %fprintf(1,'Iteration %d: %.3f\n', i, err(i));
end

end

