function [ x, res, mse, xk ] = multiGridRoughGradient( x, A, b, up, down, numIter, gt )
%MULTIGRIDROUGHGRADIENT Solve A*x = b using multigrid gradient descent.
%  x - starting point and returned solution.
%  res - If gt is given, sequence of norm(x-gt) values, else iteration
%    residuals x'(Ax-b).
%  A - square matrix.
%  b - right-hand size.
%  up - cell array of upsample operators s.t. x{i+1} = up{i}*x{i}
%  down - cell array of downsample operators s.t. x{i} = down{i}*x{i+1}
%  numIter - (optional) how many iterations to perform.
%  gt - (optional) the ground-truth solution for x, used only in
%    calculating the residual values 'res'.

numLevels = length(up)+1;
N = size(A,1);

hasxk = exist('xk','var');

res = [];
mse = [];
xk = [];

% Set default number of iterations.
if( ~exist('numIter', 'var') || isempty(numIter) )
    numIter = 100;
end

hasgt = exist('gt','var') && ~isempty(gt);

% Make the per-level A matrices.
tmp = cell(numLevels,1);
tmp{numLevels} = A;
A = tmp;
clear tmp;
F = cell(numLevels-1,1);
for i=numLevels-1:-1:1
    A{i} = down{i}*A{i+1}*up{i};
end

% REMOVE THIS!
% This replaces the initial estimate x by solving the problem directly
% on the lowest resolution, then upsampling.
%if( numLevels > 1 )
%    tmpb = b;
%    for i=numLevels-1:-1:1
%        tmpb = down{i}*tmpb;
%    end
%    x = A{1}\tmpb;
%    for i=1:numLevels-1
%        x = up{i}*x;
%    end
%end

r = cell(numLevels,1);
d = cell(numLevels,1);
k = cell(numLevels,1);
s = cell(numLevels,1);

figure;
while( true )
    numIter = numIter - 1;
    
    if( hasxk )
        xk(:, end+1) = x;
    end
    % Calculate residuals
    r{numLevels} = b - A{numLevels}*x;
    d{numLevels} = r{numLevels};
    if( hasgt )
        mse(end+1) = norm(x-gt(:),2)^2/N;
    end
    %res(end+1) = -x'*(b - 0.5*A{numLevels}*x);
    res(end+1) = norm(r{numLevels});
    
    fprintf(1, 'res: %.3f %d\n', res(end), numIter);
    if( numIter < 0 )
        break;
    end
    
    for i=numLevels-1:-1:1
        r{i} = down{i}*r{i+1};
        d{i} = r{i};
    end
    
    % Calculate level correction directions
    for i=1:numLevels
        k{i} = A{i}*d{i};
        
        for j=i:-1:2
            k{j-1} = down{j-1}*k{j};
        end
        
        s{1} = zeros(size(d{1}));
        for j=1:i-1
            s{j} = s{j} + d{j}*(k{j}'*d{j});
            s{j+1} = up{j}*s{j};
        end
        
        d{i} = d{i} - s{i};
        k{i} = A{i}*d{i};
        d{i} = d{i} / sqrt(k{i}'*d{i});
    end
    
    % Calculate correction
    s{1} = (d{1}'*r{1})*d{1};
    for i=2:numLevels
        s{i} = up{i-1}*s{i-1};
        s{i} = s{i} + d{i};
    end
    
    x = x + s{numLevels};
    %imshow(reshape(x,[618,800]));
    %drawnow;
end

end

