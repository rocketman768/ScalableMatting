function [ x, res ] = multiGridGradient3( x, A, b, up, down, numIter, gt )
%MULTIGRIDROUGHGRADIENT3 Solve A*x = b using multigrid gradient descent with rough and smooth correction directions.
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
epsilon = 1e1*eps;

res = [];

% Set default number of iterations.
if( ~exist('numIter', 'var') || isempty(numIter) )
    numIter = 100;
end

hasgt = exist('gt','var') && ~isempty(gt);

% Make the per-level A matrices and mass matrices M.
tmp = cell(numLevels,1);
tmp{numLevels} = A;
A = tmp;
clear tmp;
M = cell(numLevels,1);

M{numLevels} = speye(size(A{numLevels}));
%M{numLevels} = spdiags(spdiags(A{numLevels},0),0,size(A{numLevels}, 1),size(A{numLevels}, 2));
for i=numLevels-1:-1:1
    A{i} = down{i}*A{i+1}*up{i};
    M{i} = speye(size(A{i}));
    %M{i} = spdiags(spdiags(A{i},0),0,size(A{i}, 1),size(A{i}, 2));
end

rr = cell(numLevels,1);
rs = cell(numLevels,1);
dr = cell(numLevels,1);
ds = cell(numLevels,1);
k = cell(numLevels,1);
s = cell(numLevels,1);
g = cell(numLevels,1);
f = cell(numLevels,1);

while( true )
    numIter = numIter - 1;
    
    % Calculate rought and smooth residuals.
    rr{numLevels} = b - A{numLevels}*x;
    rs{numLevels} = M{numLevels}\rr{numLevels};
    dr{numLevels} = rr{numLevels};
    ds{numLevels} = rs{numLevels};
    if( hasgt )
        res(end+1) = norm(x-gt(:));
    else
        res(end+1) = -x'*(b - 0.5*A{numLevels}*x);
    end
    
    fprintf(1, 'res: %.3f\n', res(end));
    if( numIter < 0 )
        break;
    end
    
    for i=numLevels-1:-1:1
        rr{i} = down{i}*rr{i+1};
        rs{i} = M{i}\rr{i};
        dr{i} = rr{i};
        ds{i} = rs{i};
    end
    
    % Calculate rough and smooth level correction directions
    for i=1:numLevels
        k{i} = A{i}*dr{i};
        f{i} = A{i}*ds{i};
        
        for j=i:-1:2
            k{j-1} = down{j-1}*k{j};
            f{j-1} = down{j-1}*f{j};
        end
        
        s{1} = zeros(size(dr{1}));
        g{1} = zeros(size(dr{1}));
        for j=1:i-1
            % Eq. (2.34)
            s{j} = s{j} + dr{j}*(k{j}'*dr{j}) + ds{j}*(k{j}'*ds{j});
            % Eq. (2.35)
            g{j} = g{j} + dr{j}*(f{j}'*dr{j}) + ds{j}*(f{j}'*ds{j});
            s{j+1} = up{j}*s{j};
            g{j+1} = up{j}*g{j};
        end
        
        dr{i} = dr{i} - s{i};
        ds{i} = ds{i} - g{i};
        k{i} = A{i}*ds{i};
        % Extra sum term needs to be removed as in Eq. (2.35).
        if( i > 1 )
            s{i} = ds{i} - dr{i}*(k{i}'*dr{i});
            if( norm(s{i},Inf) > epsilon*norm(ds{i},Inf) )
                ds{i} = s{i};
            else
                ds{i} = zeros(size(s{i}));
            end
        else
            ds{i} = zeros(size(ds{i}));
        end
        
        % A-normalize rough and smooth level correction directions
        k{i} = A{i}*dr{i};
        dr{i} = dr{i}/sqrt(k{i}'*dr{i});
        if( norm(ds{i},Inf) > 0 )
            k{i} = A{i}*ds{i};
            ds{i} = ds{i}/sqrt(k{i}'*ds{i});
        end
    end
    
    % Calculate correction with rough and smooth level correction
    % directions
    s{1} = (dr{1}'*rr{1})*dr{1} + (ds{1}'*rr{1})*ds{1};
    for i=2:numLevels
        s{i} = up{i-1}*s{i-1};
        s{i} = s{i} + (dr{i}'*rr{i})*dr{i} + (ds{i}'*rr{i})*ds{i};
    end
    
    x = x + s{numLevels};
end

end

