function [ x, res, mse ] = cg( x, A, b, iter, xgt )
%CG Summary of this function goes here
%   Detailed explanation goes here

hasxgt = exist('xgt','var');

r = b-A*x;
z = r;
%z = gaussSeidel(zeros(size(x)), A, r, 3);
p = z;

res = zeros(iter,1);
if(hasxgt)
    mse = zeros(iter,1);
end

for i=1:iter
    % Restart every 50 iterations
%     if( mod(i,50) == 0 )
%         r = b-A*x;
%         z = r;
%         %z = gaussSeidel(zeros(size(x)), A, r, 3);
%         p = z;
%     end
    res(i) = norm(b-A*x);
    
    % Termination condition. Remove.
    if( res(i)/res(1) < 1e-4 )
        return;
    end
    
    if(hasxgt)
        mse(i) = norm(x-xgt(:),2)^2/size(A,1);
    end
    fprintf(1,'%.3e\n',res(i));
    alpha = r'*z/(p'*A*p);
    x = x + alpha*p;
    rnew = r - alpha*A*p;
    znew = rnew;
    %znew = gaussSeidel(zeros(size(x)), A, rnew, 3);
    beta = znew'*rnew/(z'*r);
    p = znew + beta*p;
    
    z = znew;
    r = rnew;
end

end

