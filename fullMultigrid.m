function [ x, res, mse ] = fullMultigrid( x0, A, b, up, down, iterations, xgt )
%FULLMULTIGRID Summary of this function goes here
%   Detailed explanation goes here

hasgt = exist('xgt','var');
N = size(A,1);
numLevels = length(up)+1;

Alist = cell(numLevels,1);
Alist{numLevels} = A;

for i=numLevels-1:-1:1
    Alist{i} = down{i}*Alist{i+1}*up{i};
end

x = x0;

res = zeros(iterations,1);
mse = zeros(iterations,1);
resh = figure;
for i=1:iterations
    if( hasgt )
        mse(i) = norm(x-xgt(:),2)^2/N;
    end
    res(i) = norm(b-A*x);
    
    x = fullMultigrid_internal( x, Alist, b, up, down, numLevels );
    
    fprintf(1,'%d %.1e\n',iterations-i,res(i));
    set(0,'CurrentFigure',resh);
    loglog(res);
    drawnow;
end

end

function [x] = vcycle_internal( x0, A, b, up, down, level )%+++++++++++++++

% Pre-smoothing:
x0 = gaussSeidel(x0,A{level},b,3);

% Solve on smaller grid:
r = b - A{level}*x0;
r_small = down{level-1}*r;
if level==2
%     if( size(down{1},1) < 1e4 )
%         % If the system size is small, can solve exactly in memory.
%         e = full(A{1})\r_small;
%     else
        % Otherwise, just do relaxation for a while.
        e = gaussSeidel( zeros(size(r_small)), A{1}, r_small, 10 );
        %fprintf(1,'%.2e\n',norm(A{1}*e-r_small));
%     end
else
    e = vcycle_internal( zeros(size(r_small)), A, r_small, up, down, level-1 );
end

% Fine grid correction:
x = x0 + up{level-1}*e;

% Post-smoothing:
x = gaussSeidel(x,A{level},b,1);

end%-----------------------------------------------------------------------

function [ x ] = fullMultigrid_internal( x0, A, b, up, down, level )%++++++

% I think this should be `level > 1`, but vcycle_internal below requires
% level >= 2 when called.
if( level > 2 )
    r = b - A{level}*x0;
    r_small = down{level-1}*r;
    err_small = fullMultigrid_internal(zeros(size(r_small)), A, r_small, up, down, level-1);
    err = up{level-1}*err_small;
    
    % Right?
    x = x0 + err;
    % Book says this, but I'm pretty sure it's wrong.
    %x = err;
else
    x = x0;
end

x = vcycle_internal(x, A, b, up, down, level);

end%-----------------------------------------------------------------------