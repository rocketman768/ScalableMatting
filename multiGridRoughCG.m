function [ x, res, mse, xk ] = multiGridRoughCG( x, A, b, up, down, numIter, gt, numGradIter )
%MULTIGRIDROUGHCG Solve A*x = b using multigrid conjugate gradient descent.
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

hasxk = exist('xk','var');

res = [];
mse = [];
xk = [];

% Set default number of iterations.
if( ~exist('numIter', 'var') || isempty(numIter) )
    numIter = 100;
end

hasgt = exist('gt','var') && ~isempty(gt);

% Set default gradient iterations.
if( ~exist('numGradIter', 'var') || isempty(numGradIter) )
    numGradIter = 1;
end

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
dnew = cell(numLevels,1);
k = cell(numLevels,1);
s = cell(numLevels,1);

%% First iteration is multigrid gradient descent.
for i=1:numGradIter
    numIter = numIter - 1;
    
    %xk(:,end+1) = x;
    % Calculate residuals
    r{numLevels} = b - A{numLevels}*x;
    d{numLevels} = r{numLevels};
    if( hasgt )
        % res = E[ |x-gt| ]
        %res(end+1) = norm(x-gt(:),1)/N;
        
        % MSE
        mse(end+1) = norm(x-gt(:),2)^2/N;
    end
    %res(end+1) = -x'*(b - 0.5*A{numLevels}*x);
    res(end+1) = norm(r{numLevels});
    
    
    fprintf(1, 'res: %.3f\n', res(end));
    
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
    % At this point, d{i} is \alpha_i d_i^{co} in Eq. (2.21)?
    s{1} = (d{1}'*r{1})*d{1};
    for i=2:numLevels
        s{i} = up{i-1}*s{i-1};
        s{i} = s{i} + d{i};
    end
    
    x = x + s{numLevels};
end

%% The rest is multigrid CG
%figure;
while( true )
    numIter = numIter - 1;
    
    fprintf(1, 'i: %d res: %.1e\n', numIter, res(end));
    if( numIter < 0 )
        break;
    end
    
    % Calculate residuals
    r{numLevels} = b - A{numLevels}*x;
    dnew{numLevels} = r{numLevels};
    for i=numLevels-1:-1:1
        r{i} = down{i}*dnew{i+1};
        dnew{i} = r{i};
    end
    
    if( hasxk )
        xk(:,end+1) = x;
    end
    
    if( hasgt )
        % MSE
        mse(end+1) = norm(x-gt(:),2).^2/N;
    end
    %res(end+1) = -x'*(b - 0.5*A{numLevels}*x);
    res(end+1) = norm(r{numLevels});
    
    % Termination condition. Remove.
    if( res(end)/res(1) < 1e-4 )
        return;
    end
    
    % Calculate new level correction directions
    for i=1:numLevels
        % Old direction A-orthogonal to new coarser directions
        k{i} = A{i}*d{i};
        
        for j=i:-1:2
            k{j-1} = down{j-1}*k{j};
        end
        
        s{1} = zeros(size(d{1}));
        for j=1:i-1
            s{j} = s{j} + dnew{j}*(k{j}'*dnew{j});
            s{j+1} = up{j}*s{j};
        end
        
        % A-normalize old modified level correction direction
        d{i} = d{i} - s{i};
        k{i} = A{i}*d{i};
        d{i} = d{i} / sqrt(k{i}'*d{i});
        
        % New direction A-orthogonal to new and old coarser directions
        k{i} = A{i}*dnew{i};
        for j=i:-1:2
            k{j-1} = down{j-1}*k{j};
        end
        s{1} = zeros(size(d{1}));
        for j=1:i-1
            s{j} = s{j} + dnew{j}*(k{j}'*dnew{j});
            s{j} = s{j} + d{j}*(k{j}'*d{j});
            s{j+1} = up{j}*s{j};
        end
        
        % dnew{i} now changes from w_i to \hat{d}_i^k
        dnew{i} = dnew{i} - s{i};
        
        % The following is my best guess. k{i} already equals w_i^T * A
        s{i} = dnew{i} - d{i}*(k{i}'*d{i});
        
        % New direction A-orthogonal to old modified direction
        
        % NOTE: the following two lines were in Acharya
        %k{i} = A{i}*dnew{i};
        %s{i} = dnew{i} - d{i}*(k{i}'*d{i});
        
        % NOTE: the following two lines were in Pflaum
        %k{i} = A{i}*d{i};
        %s{i} = dnew{i} - d{i}*(k{i}'*dnew{i});
        
        if( norm(s{i},Inf) > epsilon*norm(dnew{i},Inf) )
            dnew{i} = s{i};
        else
            fprintf(1,'HERE\n');
            dnew{i} = d{i};
            d{i} = zeros(size(d{i}));
        end
        
        % A-normalize new level correction direction
        k{i} = A{i}*dnew{i};
        dnew{i} = dnew{i}/sqrt(k{i}'*dnew{i});
    end
    
    % Calculate correction
    d{1} = dnew{1};
    s{1} = (d{1}'*r{1})*d{1};
    for i=2:numLevels
        d{i} = dnew{i};
        s{i} = up{i-1}*s{i-1};
        s{i} = s{i} + (d{i}'*r{i})*d{i};
    end
    
    x = x + s{numLevels};
    %imshow(reshape(x,[678,800]));
    %drawnow;
end

end
