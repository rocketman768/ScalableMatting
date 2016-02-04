% function [ x ] = vcycle( x0, A, b, up, down, level )
% %VCYCLE v-cycle recursive multigrid without smoothing
% %     x - output
% %    x0 - input
% %     A - input matrix at this level
% %    up - full list of prolongation operators
% %  down - full list of restriction operators
% % level - current level
% 
% % Pre-smoothing:
% %  x0 = x0 - W^(-1)*(A*x0-b)
% %  W := Jacobi regularizer diag(A)
% % for i=1:1
% %     x0 = x0 - (A*x0-b) ./ diag(A);
% % end
% 
% %x0 = jacobi(x0,A,b,1, 1.1);
% %fprintf(1,'Level %d>\n', level);
% x0 = gaussSeidel(x0,A,b,3);
% 
% r = b - A*x0;
% r_small = down{level}*r;
% if level==1
%     if( size(down{1},1) < 1e4 )
%         % If the system size is small, can solve exactly in memory.
%         e = full(down{1}*A*up{1})\r_small;
%     else
%         % Otherwise, just do Jacobi for a while.
%         %e = jacobi( zeros(size(r_small)), down{1}*A*up{1}, r_small, 10, 1.1 );
%         e = gaussSeidel( zeros(size(r_small)), down{1}*A*up{1}, r_small, 10 );
%     end
% else
%     e = vcycle( zeros(size(r_small)), down{level}*A*up{level}, r_small, up, down, level-1 );
% end
% 
% x = x0 + up{level}*e;
% %fprintf(1,'Level %d<\n', level);
% % Post-smoothing
% % for i=1:1
% %     x = x - (A*x-b) ./ diag(A);
% % end
% 
% %x = jacobi(x,A,b,1, 1.1);
% x = gaussSeidel(x,A,b,1);
% 
% end

function [ x, res, mse ] = vcycle( x0, A, b, up, down, iterations, xgt )%+++++++++++++++++
%VCYCLE v-cycle recursive multigrid without smoothing
%     x - output
%    x0 - input
%     A - input matrix at this level
%    up - full list of prolongation operators
%  down - full list of restriction operators
% level - current level

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
    
    % Termination condition. Remove.
    if( res(i)/res(1) < 1e-4 )
        return;
    end
    
    x = vcycle_internal( x, Alist, b, up, down, numLevels );
    
    fprintf(1,'%d %.2e\n',iterations-i,res(i));

    set(0,'CurrentFigure',resh);
    semilogy(res);
    drawnow;
end

end%-----------------------------------------------------------------------

function [x] = vcycle_internal( x0, A, b, up, down, level )%+++++++++++++++

% Pre-smoothing:
x0 = gaussSeidel(x0,A{level},b,5);

% Solve on smaller grid:
r = b - A{level}*x0;
r_small = down{level-1}*r;
if level==2
%     if( size(down{1},1) < 1e4 )
%         % If the system size is small, can solve exactly in memory.
         e = full(A{1})\r_small;
%     else
        % Otherwise, just do relaxation for a while.
        %e = gaussSeidel( zeros(size(r_small)), A{1}, r_small, 10 );
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
