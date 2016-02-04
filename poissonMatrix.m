function [ A ] = poissonMatrix( keepmask )
%POISSONMATRIX Returns a matrix for solving the poisson blending problem.
%   A*x represents x(i) if keepmask(i)==1. Otherwise, it's Laplacian(x(i)).

[m,n] = size(keepmask);
N = m*n;
idx = reshape(1:N,[m,n]);
row = repmat([1:m]',[1,n]);
col = repmat([1:n],[m,1]);

i = ones([5*N,1]);
j = ones([5*N,1]);
s = zeros([5*N,1]);

keep = idx(keepmask);
replace = idx(~keepmask);
num_replace = numel(replace);

% These parts correspond to identity matrix where we want to keep the
% original pixels.
k = numel(keep);
i(1:k) = idx(keep);
j(1:k) = idx(keep);
s(1:k) = 1;

k = k+1;

% Self.
i(k:(k+num_replace-1)) = idx(replace);
j(k:(k+num_replace-1)) = idx(replace);
s(k:(k+num_replace-1)) = 4;

k = k+num_replace;

% Top neighbor.
toprows = row(replace)-1;
centerrows = row(replace);
cols = col(replace);
valid = toprows>0;
toprows = toprows( valid );
cols = cols( valid );
centerrows = centerrows( valid );

top_indices = idx(sub2ind([m,n],toprows,cols));
center_indices = idx(sub2ind([m,n],centerrows,cols));
num_ind = numel(center_indices);

i(k:(k+num_ind-1)) = center_indices(:);
j(k:(k+num_ind-1)) = top_indices(:);
s(k:(k+num_ind-1)) = -1;

k = k+num_ind;

% Bottom neighbor.
bottomrows = row(replace)+1;
centerrows = row(replace);
cols = col(replace);
valid = bottomrows<=m;
bottomrows = bottomrows( valid );
cols = cols( valid );
centerrows = centerrows( valid );

bottom_indices = idx(sub2ind([m,n],bottomrows,cols));
center_indices = idx(sub2ind([m,n],centerrows,cols));
num_ind = numel(center_indices);

i(k:(k+num_ind-1)) = center_indices(:);
j(k:(k+num_ind-1)) = bottom_indices(:);
s(k:(k+num_ind-1)) = -1;

k = k+num_ind;

% Left neighbor.
rows = row(replace);
centercols = col(replace);
leftcols = centercols - 1;
valid = leftcols>0;
centercols = centercols( valid );
rows = rows( valid );
leftcols = leftcols( valid );

left_indices = idx(sub2ind([m,n],rows,leftcols));
center_indices = idx(sub2ind([m,n],rows,centercols));
num_ind = numel(center_indices);

i(k:(k+num_ind-1)) = center_indices(:);
j(k:(k+num_ind-1)) = left_indices(:);
s(k:(k+num_ind-1)) = -1;

k = k+num_ind;

% Right neighbor.
rows = row(replace);
centercols = col(replace);
rightcols = centercols + 1;
valid = rightcols<=n;
centercols = centercols( valid );
rows = rows( valid );
rightcols = rightcols( valid );

right_indices = idx(sub2ind([m,n],rows,rightcols));
center_indices = idx(sub2ind([m,n],rows,centercols));
num_ind = numel(center_indices);

i(k:(k+num_ind-1)) = center_indices(:);
j(k:(k+num_ind-1)) = right_indices(:);
s(k:(k+num_ind-1)) = -1;

k = k+num_ind;

A = sparse(i(1:k-1),j(1:k-1),s(1:k-1),N,N);

end
